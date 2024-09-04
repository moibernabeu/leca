#!/usr/bin/env python3

import ete3

from tree_utils import read_treeline, read_tree, read_lineage


def get_sp(header):
    return header.split('_')[0]


def main():
    ifile = 'outputs/sptree_TOLDBA/marker_trees/TOLDBA_marker_trees.nwk'
    lng = 'data/toldb.lng'
    otrees = 'outputs/marker_trees/TOLDBA_rooted_marker_trees.nwk'

    lng = read_lineage(lng)

    print(lng)

    print([lng[x]['sg'] for x in lng])

    supergroups = set([lng[x]['sg'] for x in lng])

    ostr = []
    for line in open(ifile, 'r'):
        treeline = read_treeline(line)
        tr = read_tree(treeline['tree'], lng, get_sp=get_sp,
                       root=False,
                       rooting_method='midpoint',
                       mrca_annotation=False)
        
        leaves = []
        for leaf in tr.get_leaves():
            if leaf.sg == 'Metamonada':
                leaves.append(leaf.name)

        tr.set_outgroup(leaves[0])

        for sg in set([x.sg for x in tr.get_leaves()]):
            mphy = tr.check_monophyly([sg], 'sg')
            print(treeline['tname'], sg, mphy[0], len(mphy[2]), len(tr))

    return 0


if __name__ == '__main__':
    main()