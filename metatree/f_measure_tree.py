import os
from collections import defaultdict

import dendropy
import ete3
from ete3 import Tree, TreeStyle, TextFace, add_face_to_node


class FMeasureTree(object):
    """Compares Phylorank output tables and plots differences in a tree."""

    def __init__(self, path_taxonomy):
        """Initialises the class, requires the shared taxonomy file."""
        self.files = dict()
        self.tf = TaxonomyFile(path_taxonomy)

    def add_table(self, label, path):
        """Adds a table for comparison."""
        self.files[label] = path

    def run(self, legend, out_path, rotation_deg=0):

        # Required for the script to be able to render an output.
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'

        # Create a newick tree spanning all nodes identified in f-measure tables.
        newick = NewickTree(self.tf)
        f_info = dict()
        for model_id, f_path in self.files.items():
            fm = FMeasureTable(f_path)
            newick.add_nodes(fm)
            f_info[model_id] = fm

        # Determine the data for the graph.
        f_expected = dict()
        f_exclusive_in = defaultdict(lambda: defaultdict(lambda: 0))  # [model][rank] = num
        f_exclusive_out = defaultdict(lambda: defaultdict(lambda: 0))
        f_common_in = defaultdict(lambda: 0)
        f_common_out = defaultdict(lambda: 0)

        # Index each of the ranks.
        d_rank_model = defaultdict(dict)
        for model_id, fm in f_info.items():
            for rank, rank_dict in fm.get_content().items():
                d_rank_model[rank][model_id] = rank_dict

        for rank, model_dict in d_rank_model.items():

            # This is common for all ranks.
            if set(f_info.keys()) == set(model_dict.keys()):

                # Determine the number of common taxa for all of the models.
                all_in, all_out = None, None
                for cur_model_id, cur_model_info in model_dict.items():
                    if all_in is None:
                        all_in = set(cur_model_info['g_in'])
                    else:
                        all_in = all_in.intersection(set(cur_model_info['g_in']))
                    if all_out is None:
                        all_out = set(cur_model_info['g_out'])
                    else:
                        all_out = all_out.intersection(set(cur_model_info['g_out']))

                f_common_in[rank] = len(all_in)
                f_common_out[rank] = len(all_out)

            # Determine the counts
            for cur_model_id, cur_model_info in model_dict.items():
                f_exclusive_in[cur_model_id][rank] = cur_model_info['n_in']
                f_exclusive_out[cur_model_id][rank] = cur_model_info['n_out']

                if f_expected.get(rank) != cur_model_info['n_expected']:
                    print('non expected')
                if rank not in f_expected:
                    f_expected[rank] = cur_model_info['n_expected']
                else:
                    # TODO: USe the GTDB counts
                    f_expected[rank] = min(f_expected[rank], cur_model_info['n_expected'])



        # Create an ete3 tree and annotate it
        t = Tree(str(newick), format=1, quoted_node_names=True)
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.show_scale = False
        ts.rotation = rotation_deg

        # for idx, model_id in enumerate(sorted(f_info.keys())):
        #     print(idx, model_id)

        def my_layout(node):
            MARGIN_SIZE = 3
            COLOURS = ['#ADEF29', '#F0E442', '#009E73', '#56B4E9', '#E69F00', '#911eb4', '#46f0f0', '#f032e6',
                       '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
                       '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']

            F = TextFace(node.name, tight_text=True)
            F.rotation = -rotation_deg
            F.margin_right = MARGIN_SIZE
            F.margin_left = MARGIN_SIZE
            add_face_to_node(F, node, column=0, position='branch-right')

            if node.name in ['d__Bacteria', 'd__Archaea']:

                if legend:

                    tf_expected = ete3.faces.TextFace('No. Expected', fsize=8)
                    tf_expected.background.color = '#D55E00'
                    tf_expected.margin_top = MARGIN_SIZE
                    tf_expected.margin_bottom = MARGIN_SIZE
                    tf_expected.margin_right = MARGIN_SIZE
                    tf_expected.margin_left = MARGIN_SIZE
                    tf_expected.border.width = 1
                    tf_expected.rotation = -rotation_deg
                    ete3.faces.add_face_to_node(tf_expected, node, column=2, position="branch-right")

                    tf_common_in_str = '%s\n%s' % ('No. Common In', 'No. Common Out')
                    tf_common_in = ete3.faces.TextFace(tf_common_in_str, fsize=8)
                    tf_common_in.background.color = '#CC79A7'
                    tf_common_in.margin_top = MARGIN_SIZE
                    tf_common_in.margin_bottom = MARGIN_SIZE
                    tf_common_in.margin_right = MARGIN_SIZE
                    tf_common_in.margin_left = MARGIN_SIZE
                    tf_common_in.border.width = 1
                    tf_common_in.rotation = -rotation_deg
                    ete3.faces.add_face_to_node(tf_common_in, node, column=3, position="branch-right")

                    for idx, model_id in enumerate(sorted(f_info.keys())):
                        tf_generic_str = '%s\n%s' % ('No. Unique In (%s)' % model_id, 'No. Unique Out (%s)' % model_id)
                        tf_generic = ete3.faces.TextFace(tf_generic_str, fsize=8)
                        tf_generic.background.color = COLOURS[idx]
                        tf_generic.margin_top = MARGIN_SIZE
                        tf_generic.margin_bottom = MARGIN_SIZE
                        tf_generic.margin_right = MARGIN_SIZE
                        tf_generic.margin_left = MARGIN_SIZE
                        tf_generic.border.width = 1
                        tf_generic.rotation = -rotation_deg
                        ete3.faces.add_face_to_node(tf_generic, node, column=4 + idx, position="branch-right")

            # Check if this node was not monophyletic in any of the models
            if node.name in f_expected:

                f_padd = ete3.faces.TextFace('', fsize=10)
                f_padd.margin_right = 5
                f_padd.margin_left = 5
                ete3.faces.add_face_to_node(f_padd, node, column=1, position="branch-right")

                # tf_expected = ete3.faces.TextFace(f_expected[node.name], fsize=10)
                # tf_expected.background.color = '#D55E00'
                # tf_expected.margin_top = MARGIN_SIZE
                # tf_expected.margin_bottom = MARGIN_SIZE
                # tf_expected.margin_right = MARGIN_SIZE
                # tf_expected.margin_left = MARGIN_SIZE
                # tf_expected.border.width = 1
                # tf_expected.rotation = -rotation_deg
                # ete3.faces.add_face_to_node(tf_expected, node, column=2, position="branch-right")

                tf_common_in_str = f'{f_common_in[node.name]}\n{f_common_out[node.name]}'
                tf_common_in = ete3.faces.TextFace(tf_common_in_str, fsize=10)
                tf_common_in.background.color = '#CC79A7'
                tf_common_in.margin_top = MARGIN_SIZE
                tf_common_in.margin_bottom = MARGIN_SIZE
                tf_common_in.margin_right = MARGIN_SIZE
                tf_common_in.margin_left = MARGIN_SIZE
                tf_common_in.border.width = 1
                tf_common_in.rotation = -rotation_deg
                if tf_common_in_str == '0\n0':
                    tf_common_in.opacity = 0.2
                    tf_common_in.border.width = 0
                    tf_common_in.background.color = '#f2f2f2'
                ete3.faces.add_face_to_node(tf_common_in, node, column=2, position="branch-right")

                for idx, model_id in enumerate(sorted(f_info.keys())):
                    tf_generic_str = f'{f_exclusive_in[model_id][node.name]}\n{f_exclusive_out[model_id][node.name]}\n?'
                    tf_generic = ete3.faces.TextFace(tf_generic_str, fsize=10)
                    tf_generic.background.color = COLOURS[idx]
                    tf_generic.margin_top = MARGIN_SIZE
                    tf_generic.margin_bottom = MARGIN_SIZE
                    tf_generic.margin_right = MARGIN_SIZE
                    tf_generic.margin_left = MARGIN_SIZE
                    tf_generic.border.width = 1
                    tf_generic.rotation = -rotation_deg
                    if tf_generic_str == '0\n0':
                        tf_generic.opacity = 0.2
                        tf_generic.border.width = 0
                        tf_generic.background.color = '#f2f2f2'
                    ete3.faces.add_face_to_node(tf_generic, node, column=3 + idx, position="branch-right")

        ts.layout_fn = my_layout

        # t.show(tree_style=ts)
        # xvfb-run python concat_f_measure.py

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

    def __init__(self, path):
        self.path = path
        self.content = self.read()

    def read(self):
        out = dict()
        with open(self.path) as fh:
            cols = {x: i for i, x in enumerate(fh.readline().strip().split('\t'))}
            for line in fh.readlines():
                vals = line.split('\t')

                # Only considering those which aren't monophyletic.
                f_measure = float(vals[cols['F-measure']])
                if f_measure >= 1.0:
                    continue

                taxon = vals[cols['Taxon']]
                n_expected = int(vals[cols['No. Expected in Tree']])
                rogue_in = vals[cols['Rogue in']].strip()
                rogue_out = vals[cols['Rogue out']].strip()

                hit = dict()
                hit['n_expected'] = n_expected
                hit['n_in'] = 0 if len(rogue_in) == 0 else len(rogue_in.split(','))
                hit['g_in'] = list() if len(rogue_in) == 0 else rogue_in.split(',')
                hit['n_out'] = 0 if len(rogue_out) == 0 else len(rogue_out.split(','))
                hit['g_out'] = list() if len(rogue_out) == 0 else rogue_out.split(',')
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
