import logging
import os

from metatree.common import make_sure_path_exists
from metatree.f_measure_tree import FMeasureTree
from metatree.io import Batchfile, RfResults
from metatree.io.taxonomy_file import TaxonomyFile
from metatree.tree_decorate import TreeDecorate
from metatree.tree_dist import TreeDist
from metatree.tree_root import TreeRoot


def run_pipeline(batchfile: Batchfile, out_dir: str, tax_file: TaxonomyFile, outgroup: str, cpus: int):
    logger = logging.getLogger('timestamp')

    # Setup output paths.
    dir_root = os.path.join(out_dir, 'intermediate_results', 'trees_rooted')
    dir_dec = os.path.join(out_dir, 'intermediate_results', 'trees_decorated')

    dir_rf_common = os.path.join(out_dir, 'results', 'robinson_foulds_common_taxa')
    dir_rf_all = os.path.join(out_dir, 'results', 'robinson_foulds_all_taxa')
    make_sure_path_exists(dir_rf_common)
    make_sure_path_exists(dir_rf_all)
    rf_common = RfResults(os.path.join(dir_rf_common, 'rf_common_taxa.tsv'))
    rf_all = RfResults(os.path.join(dir_rf_all, 'rf_all_taxa.tsv'))

    # tbl_diff = os.path.join(out_dir, 'results', 'model_taxonomy_diff.tsv')

    # Root the trees.
    tree_root = TreeRoot(dir_root)
    tree_root.run(batchfile, dir_root, outgroup, tax_file, cpus)

    # Decorate the trees.
    tree_decorate = TreeDecorate(dir_dec)
    tree_decorate.run(batchfile, dir_root, dir_dec, tax_file, cpus)

    # Pairwise comparison of all trees.
    td = TreeDist()
    td.run(rf_common, batchfile, dir_root, dir_dec, cpus, common_taxa=True)
    td.run(rf_all, batchfile, dir_root, dir_dec, cpus, common_taxa=False)

    # Summarise the pairwise distances (trees + heatmap)
    logger.info(f'Writing pairwise Robinson-Foulds distances for common taxa to: {dir_rf_common}')
    td.summarise_dist(rf_common, dir_rf_common)
    logger.info(f'Writing pairwise Robinson-Foulds distances for all taxa to: {dir_rf_all}')
    td.summarise_dist(rf_all, dir_rf_all)

    # Summarise the differences between all models and the reference.
    # mmt = MismatchTable(tbl_diff)
    # mmt.run_and_save(batchfile, dir_dec, tax_file)

    # Create the tree-of-trees comparison.
    fmt = FMeasureTree(tax_file.path)
    for tree_id, tree_path in batchfile.data.items():
        if tree_id != batchfile.ref:
            fmt.add_table(tree_id, os.path.join(dir_dec, f'{tree_id}_rooted_decorated.tree-table'))
    for legend in (True, False):
        if legend:
            the_path = os.path.join(out_dir, 'results', 'tree_comparison_legend.svg')
        else:
            the_path = os.path.join(out_dir, 'results', 'tree_comparison.svg')
        fmt.run(legend=legend, out_path=the_path)

    return
