import logging
import os

import dendropy

from metatree.exception import MetaTreeExit


class Batchfile(object):

    def __init__(self, path):
        self.logger = logging.getLogger('timestamp')
        self.path = path
        self.ref, self.data = self.read()

    def read(self):
        if not os.path.isfile(self.path):
            raise MetaTreeExit(f'The batchfile does not exist: {self.path}')
        out = dict()
        ref = None
        invalid_paths = list()
        with open(self.path) as fh:
            for line in fh.readlines():
                line = line.strip()
                if not line.startswith('#'):
                    tree_id, tree_path = line.split('\t')
                    out[tree_id] = tree_path
                    if ref is None:
                        ref = tree_id
                    if not os.path.isfile(tree_path):
                        invalid_paths.append((tree_id, tree_path))
        for tree_id, tree_path in invalid_paths:
            self.logger.error(f'The path for {tree_id} does not exist: {tree_path}')
        if len(invalid_paths) > 0:
            raise MetaTreeExit('Invalid tree paths were present in the batchfile.')
        return ref, out

    def common_taxa(self):
        out = set()
        for tree_id, tree_path in self.data.items():
            tree = dendropy.Tree.get_from_path(tree_path, schema='newick', preserve_underscores=True)
            cur_set = {x.label for x in tree.taxon_namespace}
            if len(out) == 0:
                out = cur_set
            else:
                out = out.intersection(cur_set)
        return out
