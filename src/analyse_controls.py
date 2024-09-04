#!/usr/bin/env python3

from tree_utils import read_lineage, read_treeline, read_tree, annotate_internal
from utils import get_list_props

import ete3


def get_sp(header):
    if 'rvdb' in header:
        return header.replace('Prok|', '').split('_', 1)[0]
    else:
        return header.split('_', 2)[1]


def analyse_sister(node: ete3.PhyloNode):
    groups = [x.sg for x in node.get_leaves()]
    props = get_list_props(groups, True)
    if (len(props) > 1 and props[0][1] == props[1][1]) or props[0][1] <= 0.5:
        return (node.mrca, 'mrca')
    return props[0]


def root_farthest(tree: ete3.PhyloTree, node: ete3.PhyloNode):
    possible_roots = []
    for st in tree.traverse():
        if (all([x not in node.get_leaf_names() for x in st.get_leaf_names()]) and
            all(['Prok|' in x for x in st.get_leaf_names()])):
            possible_roots.append((st, tree.get_distance(node, st)))
    possible_roots.sort(key=lambda x: x[1], reverse=True)

    tree.set_outgroup(possible_roots[0][0])

    return 0


def main():
    lng = read_lineage('data/LECA_proteophylum.lng')

    print('tname', 'ingroup_len', 'ingroup_mrca', 'sister_len',
          'sister_mrca', 'sister_MAB', 'sister_MAB_method',
          'ingroup_seqs_in_sister', 'ingroup_seqs', sep='\t')

    for line in open('outputs/controls/alpha_round1.nwk', 'r'):
        treeline = read_treeline(line)
        tree = read_tree(treeline['tree'], lng, root=False, get_sp=get_sp)
        mphy = tree.get_monophyletic(['Alphaproteobacteria'], 'sg')
        groups = [(group, len(group)) for group in mphy]
        
        groups.sort(key=lambda x: x[1], reverse=True)

        root_farthest(tree, groups[0][0])

        annotate_internal(tree)

        to_ingroup = True
        ingroup = groups[0][0]
        while to_ingroup and not ingroup.up.is_root():
            sister = ingroup.get_sisters()[0]
            prop_ogseqs = sum([True for x in sister.get_leaf_names() if 'Prok|' not in x]) / len(sister)
            if prop_ogseqs > 0.3:
                to_ingroup = True
                ingroup = ingroup.up
            else:
                to_ingroup = False
        
        sister = ingroup.get_sisters()[0]
        sister_analysis = analyse_sister(sister)

        prop_ogseqs = sum([True for x in sister.get_leaf_names() if 'Prok|' not in x]) / len(sister)

        ingroup_seqs = [x for x in ingroup.get_leaf_names() if 'Prok|' not in x]

        print(treeline['tname'], len(ingroup), ingroup.mrca, len(sister),
              sister.mrca, sister_analysis[0], sister_analysis[1],
              prop_ogseqs, ';'.join(ingroup_seqs), sep='\t')

    return 0


if __name__ == '__main__':
    main()