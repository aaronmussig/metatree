import logging
import os
import tempfile
from collections import defaultdict
from multiprocessing import Pool
from warnings import simplefilter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.cluster.hierarchy import ClusterWarning
from tqdm import tqdm

from metatree.external.tree_compare import TreeCompare
from metatree.io import Batchfile, RfResults

simplefilter("ignore", ClusterWarning)


class TreeDist(object):

    def __init__(self):
        self.logger = logging.getLogger('timestamp')

    @staticmethod
    def worker(task):
        tid_a, path_a, tid_b, path_b, set_common = task

        if set_common:
            with tempfile.TemporaryDirectory() as tmp_dir:
                taxa_list = os.path.join(tmp_dir, 'taxa_list.tsv')
                with open(taxa_list, 'w') as fh:
                    for taxa in set_common:
                        fh.write(f'{taxa}\n')

                # Calculate the RF distance.
                tc = TreeCompare()
                rf, norm_rf = tc.robinson_foulds(path_a, path_b, taxa_list)

        else:
            # Calculate the RF distance.
            tc = TreeCompare()
            rf, norm_rf = tc.robinson_foulds(path_a, path_b, None)

        return tid_a, tid_b, rf, norm_rf

    def run(self, rf_results: RfResults, batchfile: Batchfile, dir_root, dir_dec, cpus: int, common_taxa: bool):

        # Determine if a common subset of taxa should be used.
        if common_taxa:
            set_common = batchfile.common_taxa()
            self.logger.info(f'Robinson-Foulds metrics will only consider those {len(set_common):,} '
                             f'taxa which are common between ALL trees.')
        else:
            set_common = None

        # Determine which trees still need to be processed.
        queue = list()
        tree_ids = list(batchfile.data.keys())
        for i in range(len(tree_ids)):
            tid_a = tree_ids[i]
            path_a = batchfile.data[tid_a]
            for j in range(i):
                tid_b = tree_ids[j]
                path_b = batchfile.data[tid_b]
                if not rf_results.is_done(tid_a, tid_b):
                    queue.append((tid_a, path_a, tid_b, path_b, set_common))

        self.logger.info(f'Calculating Robinson-Foulds distances.')
        with Pool(processes=cpus) as pool:
            for tid_a, tid_b, rf, norm_rf in tqdm(pool.imap_unordered(TreeDist.worker, queue), total=len(queue)):
                rf_results.add(tid_a, tid_b, rf, norm_rf)

        rf_results.write()

    def summarise_dist(self, rf_results: RfResults, dir_out):

        for use_norm in (True, False):
            if use_norm:
                path_out = os.path.join(dir_out, 'rf_normed.tree')
                path_hm = os.path.join(dir_out, 'rf_normed_heatmap.svg')
                plt_title = 'Normalised Robinson-Foulds Distance'
            else:
                path_out = os.path.join(dir_out, 'rf_un_normed.tree')
                path_hm = os.path.join(dir_out, 'rf_un_normed_heatmap.svg')
                plt_title = '(un)Normalised Robinson-Foulds Distance'

            metrics = defaultdict(dict)
            names = set()
            for (tid_a, tid_b), (rf, norm_rf) in rf_results.data.items():
                if use_norm:
                    metrics[tid_a][tid_b] = norm_rf
                    metrics[tid_b][tid_a] = norm_rf
                else:
                    metrics[tid_a][tid_b] = rf
                    metrics[tid_b][tid_a] = rf
                names.add(tid_a)
                names.add(tid_b)

            labels = sorted(list(names))
            mat_vals = list()
            mat = np.zeros((len(labels), len(labels)))
            for i in range(len(labels)):
                cur_row = list()
                tid_a = labels[i]
                for j in range(i + 1):
                    tid_b = labels[j]
                    if tid_a == tid_b:
                        cur_row.append(0.0)
                    else:
                        cur_row.append(metrics[tid_a][tid_b])
                        mat[i, j] = metrics[tid_a][tid_b]
                mat_vals.append(cur_row)
            mat = mat + mat.T

            # Newick
            dm = DistanceMatrix(names=labels, matrix=mat_vals)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)

            Phylo.write(tree, path_out, 'newick')

            # Heatmap
            cmap = sns.cubehelix_palette(100, reverse=True)

            sns.set(font_scale=1)
            fig_size = (15, 15)

            rf_df = pd.DataFrame(mat, columns=labels, index=labels)
            sns.clustermap(rf_df, annot=True, fmt='.3f', cmap=cmap, figsize=fig_size).fig.suptitle(plt_title)
            plt.savefig(path_hm)
