#!/usr/bin/env python3

import ete3
from LECA_utils import read_lineage, get_sp, is_LECA
from LECA_family_origin import get_farthest_neuk_node, get_sisters, get_leca, get_list_props
from tree_utils import read_treeline, read_tree
import pandas as pd
import argparse as ap

def analyse_viruses(tr, node):
    sisters = get_sisters(tr, node)
    
    if len(sisters) > 0:
        odict = {}
        odict['ntips'] = len(tr)
        odict['leca_tips'] = len(node)
        odict['sister_tips'] = len(tr) - len(node)
        groups = get_list_props([(x.d, x.sg) for x in sisters[0].get_leaves()], True)

        if groups[0][0][0] == 'Viruses':
            viruses = set([groups[0][0][1]])
            i = 0
            nvirseqs = 0
            while groups[0][0][0] == 'Viruses' and i + 1 < len(sisters):
                viruses.add(groups[0][0][1])
                nvirseqs += len(sisters[i])
                i += 1
                groups = get_list_props([(x.d, x.sg) for x in sisters[i].get_leaves()], True)

            odict['viruses'] = ';'.join(viruses)
            odict['nvirseqs'] = nvirseqs
            if groups[0][0][0] != 'Viruses' and i > 0:
                odict['nv_sister_domain'] = groups[0][0][0]
                odict['nv_sister'] = groups[0][0][1]
                odict['nv_sister_size'] = len(sisters[i])

            return odict
    return None


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-t', '--trees', dest='treefile',
                        help='Treefile.')
    parser.add_argument('-l', '--lineage', dest='lineage',
                        help='Lineage file.')
    parser.add_argument('-o', '--output', dest='ofile',
                        help='Output file in .tsv format')
    args = parser.parse_args()

    treefile = args.treefile
    lineage = args.lineage
    ofile = args.ofile

    lng = read_lineage(lineage)

    odfl = list()
    for line in open(treefile, 'r'):
        treeline = read_treeline(line)
        tr = read_tree(treeline['tree'], lng, root=False, get_sp=get_sp)

        leca = get_leca(tr)

        i = 1
        while leca is not None:
            to_keep = [x for x in tr.get_leaf_names() if x not in leca[1].get_leaf_names()]
            if leca[0]:
                leca_analysis = analyse_viruses(tr, leca[1])
                if leca_analysis is not None:
                    criteria = is_LECA([x.name for x in leca[1].get_leaves() if x.d == 'Eukaryota'])
                    odfl.append({'og': treeline['tname'], 'LECA': f'LECA{i}', **criteria, **leca_analysis})
            tr.prune(to_keep)
            leca = get_leca(tr)
        print(treeline['tname'])
    odf = pd.DataFrame(odfl)
    odf.to_csv(ofile, sep='\t', index=False)
    return 0


if __name__ == '__main__':
    main()