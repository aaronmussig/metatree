import logging
import os
import subprocess
from multiprocessing import Pool

from genometreetk import __version__ as genometreetk_v
from tqdm import tqdm

from metatree.common import make_sure_path_exists
from metatree.exception import MetaTreeExit
from metatree.io import Batchfile
from metatree.io.taxonomy_file import TaxonomyFile


class TreeRoot(object):

    def __init__(self, dir_root):
        self.dir_root = dir_root
        self.logger = logging.getLogger('timestamp')
        make_sure_path_exists(dir_root)

    @staticmethod
    def worker(task):
        tree_id, tree_in, tree_out, outgroup, tax_file = task
        args = ['genometreetk', 'outgroup', tree_in, tax_file, outgroup, tree_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        proc.communicate()
        if not proc.returncode == 0:
            raise MetaTreeExit(f'Non-zero return code: {" ".join(args)}')

    def run(self, batchfile: Batchfile, dir_root: str, outgroup: str, tax_file: TaxonomyFile, cpus: int):
        queue = list()
        for tree_id, tree_in in batchfile.data.items():
            tree_out = os.path.join(dir_root, f'{tree_id}_rooted.tree')
            if not os.path.isfile(tree_out):
                queue.append((tree_id, tree_in, tree_out, outgroup, tax_file.path))

        self.logger.info(f'Rooting trees using GenomeTreeTk v{genometreetk_v}')
        with Pool(processes=cpus) as pool:
            for _ in tqdm(pool.imap_unordered(TreeRoot.worker, queue), total=len(queue)):
                pass
