import logging
import os

from metatree.common import make_sure_path_exists
from metatree.exception import MetaTreeExit


class RfResults(object):

    def __init__(self, path):
        self.logger = logging.getLogger('timestamp')
        self.path = path
        self.data = self.read()
        make_sure_path_exists(os.path.dirname(path))

    def read(self):
        out = dict()
        if os.path.isfile(self.path):
            with open(self.path) as fh:
                for line in fh.readlines():
                    tid_a, tid_b, rf, norm_rf = line.strip().split('\t')
                    out[(tid_a, tid_b)] = (float(rf), float(norm_rf))
        return out

    def is_done(self, tid_a, tid_b):
        return (tid_a, tid_b) in self.data or (tid_b, tid_a) in self.data

    def add(self, tid_a, tid_b, rf, norm_rf):
        self.data[(tid_a, tid_b)] = (rf, norm_rf)

    def write(self):
        done = dict()
        with open(self.path, 'w') as fh:
            for (tid_a, tid_b), (rf, norm_rf) in self.data.items():
                if (tid_a, tid_b) in done and self.data[(tid_a, tid_b)] != done[(tid_a, tid_b)]:
                    raise MetaTreeExit('Inconsistent results, report this issue.')
                if (tid_b, tid_a) in done and self.data[(tid_b, tid_a)] != done[(tid_b, tid_a)]:
                    raise MetaTreeExit('Inconsistent results, report this issue.')
                if (tid_a, tid_b) in done or (tid_b, tid_a) in done:
                    continue
                else:
                    fh.write(f'{tid_a}\t{tid_b}\t{rf}\t{norm_rf}\n')
                    done[(tid_a, tid_b)] = (rf, norm_rf)
        self.logger.info(f'Pairwise Robinson-Foulds distances written to: {self.path}')
