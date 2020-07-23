import logging
import os

from metatree.exception import MetaTreeExit


class TaxonomyFile(object):

    def __init__(self, path):
        self.path = path
        self.data = self.read()
        self.logger = logging.getLogger('timestamp')

    def read(self):
        if not os.path.isfile(self.path):
            raise MetaTreeExit(f'The taxonomy file does not exist: {self.path}')
        out = dict()
        invalid_tax = list()
        with open(self.path) as fh:
            for line in fh.readlines():
                gid, tax = line.strip().split('\t')
                out[gid] = tax
                if len(tax.split(';')) != 7:
                    invalid_tax.append((gid, tax))
        for gid, tax in invalid_tax:
            self.logger.error(f'Genome {gid} does not have a valid taxonomy: {tax}')
        if len(invalid_tax) > 0:
            raise MetaTreeExit(f'There were {len(invalid_tax)} invalid taxonomies detected.')
        return out
