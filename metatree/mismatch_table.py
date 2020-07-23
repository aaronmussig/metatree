import logging
import os
from collections import defaultdict

import pandas as pd

from metatree.io import Batchfile


class FMeasureTable(object):

    def __init__(self, path):
        self.path = path
        self.content = self.read()

    def read(self):
        out = dict()

        df = pd.read_csv(self.path, delimiter='\t')
        df_filtered = df[df['F-measure'] < 1.0]

        for idx, row in df_filtered.iterrows():
            hit = dict()
            hit['n_expected'] = row['No. Expected in Tree']
            hit['n_in'] = 0 if str(row['Rogue in']) == 'nan' else len(row['Rogue in'].split(','))
            hit['g_in'] = list() if str(row['Rogue in']) == 'nan' else row['Rogue in'].split(',')
            hit['n_out'] = 0 if str(row['Rogue out']) == 'nan' else len(row['Rogue out'].split(','))
            hit['g_out'] = list() if str(row['Rogue out']) == 'nan' else row['Rogue out'].split(',')
            out[row['Taxon']] = hit

        return out

    def get_content(self) -> dict:
        return self.content


def parse_tax_file(path):
    out = dict()
    with open(path, 'r') as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax
    return out


def tax_to_dict(tax):
    out = defaultdict(list)
    for rank in tax.split(';'):
        out[rank[0]].append(rank)
    return out


def compare_tax(truth, model, comp):
    n_disagree = 0
    n_poly = 0
    n_agree = 0

    for gid, m_tax in comp.items():
        t_tax = truth[gid]

        t_tax, m_tax = t_tax.replace('; ', ';'), m_tax.replace('; ', ';')
        dict_t_tax, dict_m_tax = tax_to_dict(t_tax), tax_to_dict(m_tax)

        # Case 1: Agreement
        if m_tax == t_tax:
            n_agree += 1
            continue

        for rank in ['d', 'p', 'c', 'o', 'f', 'g', 's']:
            m_rank, t_rank = dict_m_tax[rank], dict_t_tax[rank]

            if m_rank == t_rank:
                continue

            # Case 2: Disagreement
            elif len(m_rank) == len(t_rank):
                n_disagree += 1

                print(t_tax)
                print(m_tax)
                print('--')

                break

            # Case 3: Polyphyletic / novel rank
            else:
                n_poly += 1
                break

    print(f'{model}\t{n_agree}\t{n_disagree}\t{n_poly}')
    return model, n_agree, n_disagree, n_poly


class MismatchTable(object):

    def __init__(self, path):
        self.path = path
        self.logger = logging.getLogger('timestamp')

    def run_and_save(self, batchfile: Batchfile, dir_decorated, path_tax):
        truth_tax = parse_tax_file(path_tax.path)

        for tree_id, tree_path in batchfile.data.items():
            cur_tax = parse_tax_file(os.path.join(dir_decorated, f'{tree_id}_rooted_decorated.tree-taxonomy'))
            compare_tax(truth_tax, tree_id, cur_tax)

        return
