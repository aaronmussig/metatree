import os
from collections import defaultdict

import dendropy
import ete3

class FMeasureTree(object):
    """Compares Phylorank output tables and plots differences in a tree."""

    def __init__(self, path_taxonomy):
        """Initialises the class, requires the shared taxonomy file."""
        self.files = dict()
        self.tf = TaxonomyFile(path_taxonomy)

    def add_table(self, label, path):
        """Adds a table for comparison."""
        self.files[label] = FMeasureTable(path)

    def get_n_common(self):
        """Determine the number of genomes which are common between ALL models."""
        out = defaultdict(lambda: defaultdict(lambda: 0))

        # Index by rank.
        d_rank_model = defaultdict(dict)
        for model_id, fm in self.files.items():
            for rank, rank_dict in fm.get_content().items():
                d_rank_model[rank][model_id] = rank_dict

        # Iterate over each rank and find the number of common taxa.
        for rank, model_dict in d_rank_model.items():

            # This is common for all ranks.
            if set(self.files.keys()) == set(model_dict.keys()):

                # Determine the number of common taxa for all of the models.
                all_in, all_out = None, None
                for cur_model_id, cur_model_info in model_dict.items():
                    if all_in is None:
                        all_in = set(cur_model_info['rogue_in'])
                    else:
                        all_in = all_in.intersection(set(cur_model_info['rogue_in']))
                    if all_out is None:
                        all_out = set(cur_model_info['rogue_out'])
                    else:
                        all_out = all_out.intersection(set(cur_model_info['rogue_out']))

                out[rank]['in'] = len(all_in)
                out[rank]['out'] = len(all_out)
        return out

    def get_poly_ranks(self):
        out = set()
        for model_id, fm in self.files.items():
            for rank, rank_dict in fm.get_content().items():
                if rank_dict['f_measure'] < 1.0:
                    out.add(rank)
        return out

    @staticmethod
    def get_text_face(text, bg, opacity, fsize=8, bold=False):
        margin_size = 2

        tf = ete3.faces.TextFace(str(text), fsize=fsize, tight_text=True)
        tf.background.color = bg
        tf.margin_top = margin_size
        tf.margin_bottom = margin_size
        tf.margin_right = margin_size
        tf.margin_left = margin_size
        tf.opacity = opacity
        tf.rotation = -0
        tf.bold = bold
        tf.border.color = '#f2f2f2'
        tf.border.width = None  # try 1?
        tf.hz_align = 1
        tf.vt_align = 1
        return tf

    def run(self, legend, out_path, rotation_deg=0):

        # Required for the script to be able to render an output.
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'

        # Create a newick tree spanning all nodes identified in f-measure tables.
        newick = NewickTree(self.tf)
        [newick.add_nodes(fm) for fm in self.files.values()]

        # Determine the counts for the common taxa.
        d_rank_common = self.get_n_common()

        # Create an ete3 tree and annotate it
        t = ete3.Tree(str(newick), format=1, quoted_node_names=True)
        ts = ete3.TreeStyle()
        ts.show_leaf_name = False
        ts.show_scale = False
        ts.rotation = rotation_deg

        def my_layout(node):
            COLOURS = ['#ADEF29', '#F0E442', '#009E73', '#56B4E9', '#E69F00', '#911eb4', '#46f0f0', '#f032e6',
                       '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
                       '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']

            F = ete3.TextFace(node.name, fsize=8, tight_text=True)
            F.rotation = -rotation_deg
            F.margin_right = 2
            F.margin_left = 3
            ete3.add_face_to_node(F, node, column=0, position='branch-right')

            tf_invis = FMeasureTree.get_text_face('-', None, 0.0)

            if node.name in ['d__Bacteria', 'd__Archaea']:

                # Attach a legend to the plot.
                if legend:
                    tf_legend_common_in = FMeasureTree.get_text_face('No. Common In', '#F2F2F2', 0.5, fsize=7,
                                                                     bold=True)
                    ete3.faces.add_face_to_node(tf_legend_common_in, node, column=1, position="branch-top")

                    tf_legend_common_out = FMeasureTree.get_text_face('No. Common Out', '#F2F2F2', 0.5, fsize=7,
                                                                      bold=True)
                    ete3.faces.add_face_to_node(tf_legend_common_out, node, column=1, position="branch-top")

                    tf_legend_common_blank = FMeasureTree.get_text_face('X', None, 0.0, fsize=7)
                    ete3.faces.add_face_to_node(tf_legend_common_blank, node, column=1, position="branch-top")

                    for idx, model_id in enumerate(sorted(self.files.keys())):
                        tf_legend_in = FMeasureTree.get_text_face(f'No. Rogue In ({model_id})', COLOURS[idx], 1.0,
                                                                  fsize=7)
                        ete3.faces.add_face_to_node(tf_legend_in, node, column=2 + idx, position="branch-top")

                        tf_legend_out = FMeasureTree.get_text_face(f'No. Rogue Out ({model_id})', COLOURS[idx], 1.0,
                                                                   fsize=7)
                        ete3.faces.add_face_to_node(tf_legend_out, node, column=2 + idx, position="branch-top")

                        tf_legend_exp = FMeasureTree.get_text_face(f'No. Expected ({model_id})', "#F2F2F2", 0.2,
                                                                   fsize=7)
                        ete3.faces.add_face_to_node(tf_legend_exp, node, column=2 + idx, position="branch-top")

            # Check if this node was not monophyletic in any of the models
            poly_ranks = self.get_poly_ranks()
            if node.name in poly_ranks:

                # Create a spacer between the rank name and the values.
                ete3.faces.add_face_to_node(tf_invis, node, column=1, position="branch-right")

                n_common_in = d_rank_common[node.name]['in']
                n_common_out = d_rank_common[node.name]['out']

                # Add a spacer to act as a top margin.
                ete3.faces.add_face_to_node(tf_invis, node, column=2, position="branch-right")

                # Add the number of common taxa in.
                if n_common_in == 0:
                    tf_common_in = FMeasureTree.get_text_face(n_common_in, '#f2f2f2', 0.2, bold=True)
                else:
                    tf_common_in = FMeasureTree.get_text_face(n_common_in, '#f2f2f2', 0.5, bold=True)
                ete3.faces.add_face_to_node(tf_common_in, node, column=2, position="branch-right")

                # Add the number of common taxa out.
                if n_common_out == 0:
                    tf_common_out = FMeasureTree.get_text_face(n_common_out, '#f2f2f2', 0.2, bold=True)
                else:
                    tf_common_out = FMeasureTree.get_text_face(n_common_out, '#f2f2f2', 0.5, bold=True)
                ete3.faces.add_face_to_node(tf_common_out, node, column=2, position="branch-right")

                tf_common_blank = FMeasureTree.get_text_face('X', None, 0.0)
                ete3.faces.add_face_to_node(tf_common_blank, node, column=2, position="branch-right")

                # Add model specific information for this current rank.
                for idx, model_id in enumerate(sorted(self.files.keys())):
                    cur_f_table = self.files[model_id].content
                    n_rogue_in = len(cur_f_table[node.name]["rogue_in"])
                    n_rogue_out = len(cur_f_table[node.name]["rogue_out"])

                    # Add a top margin spacer.
                    ete3.faces.add_face_to_node(tf_invis, node, column=3 + idx, position="branch-right")

                    # Number of rogue taxa in.
                    if n_rogue_in == 0:
                        tf_generic_in = FMeasureTree.get_text_face(n_rogue_in, '#f2f2f2', 0.2)
                    else:
                        tf_generic_in = FMeasureTree.get_text_face(n_rogue_in, COLOURS[idx], 1.0)
                    ete3.faces.add_face_to_node(tf_generic_in, node, column=3 + idx, position="branch-right")

                    # Number of rogue taxa out.
                    if n_rogue_out == 0:
                        tf_generic_out = FMeasureTree.get_text_face(n_rogue_out, '#f2f2f2', 0.2)
                    else:
                        tf_generic_out = FMeasureTree.get_text_face(n_rogue_out, COLOURS[idx], 1.0)
                    ete3.faces.add_face_to_node(tf_generic_out, node, column=3 + idx, position="branch-right")

                    # Expected number of taxa in this model.
                    tf_expected_cnt = FMeasureTree.get_text_face(cur_f_table[node.name]["n_expected"], '#f2f2f2', 0.2)
                    ete3.faces.add_face_to_node(tf_expected_cnt, node, column=3 + idx, position="branch-right")

        ts.layout_fn = my_layout

        if out_path is None:
            print('Run as: xvfb-run python concat_f_measure.py')
            t.show(tree_style=ts)
        else:
            t.render(out_path, tree_style=ts)


class TaxonomyFile(object):

    def __init__(self, path):
        self.path = path
        self.contents = self.read()

    def read(self):
        out = dict()
        with open(self.path, 'r') as f:
            for line in f.readlines():
                cols = line.strip().split('\t')
                accession = cols[0]
                ranks = cols[1].split(';')
                out[accession] = ranks
        return out

    def get_content(self):
        return self.contents

    def get_taxon_namespace(self):
        out = set()

        for v in self.contents.values():
            [out.add(x) for x in v]

        return out

    def get_ranks_above(self, search_str):
        for v in self.contents.values():
            last_ranks = list()
            for rank in v:
                last_ranks.append(rank)
                if rank == search_str:
                    return last_ranks


class Node(object):

    def __init__(self, name, parent):
        self.name = name
        self.children = list()
        self.parent = parent
        self.count = None
        self.metadata = list()

    def get_name(self):
        return self.name

    def add_child(self, node):
        self.children.append(node)

    def get_child(self, name):
        for child in self.children:
            if child.get_name() == name:
                return child
        return None

    def add_ranks(self, ranks):
        last_node = self
        for rank in ranks:
            search_node = last_node.get_child(rank)
            if not search_node:
                new_node = Node(rank, parent=last_node)
                last_node.add_child(new_node)
                last_node = new_node
            else:
                last_node = search_node
        return last_node

    def add_expected_count(self, count):
        self.count = count

    def __str__(self):

        parent_str = 'null' if self.parent is None else self.parent.get_name()

        expected_str = '' if self.count is None else ',\n"n_expected": %s' % self.count
        metadata_str = '' if len(self.metadata) == 0 else ',\n"metadata": [\n%s]' % '\n'.join(self.metadata)

        children_str = str()
        if len(self.children) > 0:
            children_str = ''',
            "children": [
            %s
            ]''' % ', '.join([str(x) for x in self.children])

        out = '''{
        "name": "%s",
        "parent": "%s" %s %s %s 
    }
        ''' % (self.name, parent_str, children_str, expected_str, metadata_str)

        return out

    def add_metadata(self, item):
        self.metadata.append(item)


class FMeasureTable(object):
    cols = ('Taxon', 'No. Expected in Tree', 'F-measure', 'Precision', 'Recall',
            'No. Genomes from Taxon', 'No. Genome In Lineage', 'Rogue out',
            'Rogue in')

    def __init__(self, path):
        self.path = path
        self.content = self.read()

    def read(self):
        out = dict()
        with open(self.path) as fh:
            read_cols = tuple([x for x in fh.readline().strip().split('\t')])
            if self.cols != read_cols:
                raise Exception('PhyloRank output file has different headers.')
            col_ids = {x: i for i, x in enumerate(self.cols)}

            for line in fh.readlines():
                vals = [x.strip() for x in line.split('\t')]

                taxon = vals[col_ids['Taxon']]
                n_expected = int(vals[col_ids['No. Expected in Tree']])
                f_measure = float(vals[col_ids['F-measure']])
                precision = float(vals[col_ids['Precision']])
                recall = float(float(vals[col_ids['Recall']]))
                n_from_taxon = int(vals[col_ids['No. Genomes from Taxon']])
                n_from_lineage = int(vals[col_ids['No. Genome In Lineage']])
                rogue_in = vals[col_ids['Rogue in']]
                rogue_out = vals[col_ids['Rogue out']]

                hit = dict()
                hit['n_expected'] = n_expected
                hit['f_measure'] = f_measure
                hit['precision'] = precision
                hit['recall'] = recall
                hit['n_from_taxon'] = n_from_taxon
                hit['n_from_lineage'] = n_from_lineage
                hit['rogue_in'] = list() if len(rogue_in) == 0 else rogue_in.split(',')
                hit['rogue_out'] = list() if len(rogue_out) == 0 else rogue_out.split(',')
                out[taxon] = hit
        return out

    def get_content(self) -> dict:
        return self.content


class NewickTree(object):

    def __init__(self, tf: TaxonomyFile):
        self.tf = tf
        self.taxon_namespace = dendropy.TaxonNamespace(self.tf.get_taxon_namespace())
        self.tree = dendropy.Tree()
        self.root = self.tree.seed_node

    def add_nodes(self, fm: FMeasureTable):
        for taxon, f_dict in fm.get_content().items():

            # Only constructing nodes which are polyphyletic.
            if f_dict['f_measure'] >= 1.0:
                continue

            last_node = self.root
            for cur_rank in self.tf.get_ranks_above(taxon):
                cur_node = self.tree.find_node_with_taxon_label(cur_rank)
                if not cur_node:
                    if cur_rank in ['d__Archaea', 'd__Bacteria']:
                        self.root.taxon = self.taxon_namespace.get_taxon(cur_rank)
                        cur_node = self.root
                    else:
                        new_node = dendropy.Node()
                        last_node.add_child(new_node)
                        new_node.taxon = self.taxon_namespace.get_taxon(cur_rank)
                        cur_node = new_node
                last_node = cur_node

    def __str__(self):
        return self.tree.as_string('newick')
