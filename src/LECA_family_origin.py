#!/usr/bin/env python3

from tree_utils import read_lineage, read_treeline, read_tree
from LECA_utils import get_sp, is_LECA
from utils import get_list_props, create_folder
import statistics
import pandas as pd
import argparse as ap

import ete3
import sys


def get_farthest_neuk_node(tree: ete3.PhyloTree, node: ete3.PhyloNode=None):
    if node is None:
        node = tree.get_tree_root()

    roots = []
    not_leaves = 0
    for st in tree.traverse():
        if (all([x.d != 'Eukaryota' for x in st.get_leaves()])):
            roots.append((tree.get_distance(st, node, topology_only=True), st.is_leaf(), st))
            if not st.is_leaf():
                not_leaves += 1

    roots.sort(key=lambda x: x[0], reverse=True)
    if not_leaves > 0:
        roots = [x for x in roots if not x[1]]

    return roots[0][2]


def get_sisters(tree: ete3.PhyloTree, node: ete3.PhyloNode):
    if node not in tree:
        return None

    sisters = []
    sisters.append(node.get_sisters()[0])
    upnode = sisters[0].up
    while not upnode.is_root():
        sisters.append(upnode.get_sisters()[0])
        upnode = upnode.up

    return sisters


def get_leca(tree: ete3.PhyloTree):
    # Getting Eukaryotic monophyletic groups
    emps = list(tree.get_monophyletic(values=['Eukaryota'], target_attr='d'))

    if len(emps) == 0:
        return None

    # Getting the biggest eukaryotic monophyletic group
    biggest = [(len(x), x) for x in emps]
    biggest.sort(key = lambda x: x[0], reverse=True)
    biggest = biggest[0][1]

    # Setting the outgroup in the the farthest non-eukaryotic internal node
    og = get_farthest_neuk_node(tree, biggest)
    tree.set_outgroup(og)

    # Getting the sisters
    sisters = get_sisters(tree, biggest)
    sist_summ = []
    for i, sister in enumerate(sisters):
        props = get_list_props([x.d for x in sister.get_leaves()])
        neuk_sgr_count = len(set([x.g for x in sister.get_leaves() if x.d != 'Eukaryota']))
        sist_size = len(sister)
        if 'Eukaryota' in props.keys():
            neuk_prop = 1 - props['Eukaryota']
        else:
            neuk_prop = 1

        sist_summ.append({**props, 'neuk_sgr_count': neuk_sgr_count,
                          'neuk_prop': neuk_prop, 'size': sist_size})

    # Iterate through sisters incorporating them just if they fulfil
    # the conditions
    in_leca = True
    leca = biggest
    i = 0
    while in_leca and not leca.up.is_root():
        s1 = sist_summ[i]
        try:
            s2 = sist_summ[i + 1]
        except IndexError:
            s2 = {'neuk_prop': 1}

        # The first sister is mostly eukaryotic, the threshold is stronger for
        # bigger groups
        condition_s1 = ((s1['neuk_prop'] <= 0.5 and s1['size'] <= 20) or
                        (s1['neuk_prop'] <= 0.25 and s1['size'] > 20))
        
        # The second sister is eukaryotic when the first sister is not bigger
        # than 75% of the size of the current LECA node and the first sister
        # has no more than 2 non-eukaryotic phyla in it
        condition_s2 = (s2['neuk_prop'] <= 0.25 and
                        # s2['size'] > 0.05 * len(leca) and # Excluded condition
                        s1['size'] <= 0.75 * len(leca) and
                        s1['neuk_sgr_count'] < 3)

        if (condition_s1 or condition_s2):
            leca = leca.up
        else:
            in_leca = False

        i += 1

    leca_euks = [x.name for x in leca.get_leaves() if x.d == 'Eukaryota']
    criteria = is_LECA(leca_euks)

    # Annotate as a LECA group and return it if it fulfils the two stems
    # criteria and three or five supergroups
    if criteria['is_LECA']:
        leca.add_feature('name', 'LECA')
        return [True, leca]
    elif len(leca) > 0:
        return [False, leca]
    else:
        return None


def analyse_sisters(tr: ete3.PhyloTree, node: ete3.PhyloNode,
                    group2dom: dict, max_sists=None):
    if max_sists is None:
        max_sists = 10

    sisters = get_sisters(tr, node)
    odict = {}
    for i, sister in enumerate(sisters):
        if len(odict) <= max_sists:
            domains = [x.d for x in sister.get_leaves()]

            domains = get_list_props(domains, sort=True)
            if 'Eukaryota' in dict(domains).keys():
                if dict(domains)['Eukaryota'] > 0.25:
                    return None

            groups = [x.sg for x in sister.get_leaves() if x.d != 'Eukaryota']
            groups = get_list_props(groups, sort=True)

            try:
                mixed = groups[0][1] == groups[1][1]
            except IndexError:
                mixed = False

            try:
                s2_groups = [x.g for x in sisters[i + 1].get_leaves()]
                s1_in_s2 = groups[0][0] in s2_groups
            except IndexError:
                s1_in_s2 = False

            if mixed:
                if len(set(dict(groups).keys())) > 1:
                    doms = sorted([x[0] for x in domains if x[0] != 'Eukaryota'])
                    if 'Proteobacteria' in doms and 'Bacteria' in doms:
                        doms.remove('Proteobacteria')
                    mrca = ';'.join(doms)
                    MAB_group_domain = mrca
                else:
                    mrca = groups[0][0]
                    MAB_group_domain = group2dom[groups[0][0]]
            else:
                mrca = groups[0][0]
                MAB_group_domain = group2dom[groups[0][0]]

            sdict = {'size': len(sister),
                     'groups': dict(groups),
                     'domains': dict(domains),
                     'MAB_group': groups[0][0],
                     'mrca': mrca,
                     'MAB_group_domain': MAB_group_domain,
                     'MAB_domain': domains[0][0],
                     'mixed': mixed,
                     'support': sister.support,
                     'seqs': ';'.join(sister.get_leaf_names()),
                     f'S{i + 1}_in_S{i + 2}': s1_in_s2}

            odict[f'S{i + 1}'] = sdict
        else:
            break

    return odict


def prune_longbranchozoa(node: ete3.PhyloNode, maxdev=2):
    node2tip = {x.name: node.get_distance(x) for x in node.get_leaves() if x.d == 'Eukaryota'}
    median = statistics.median(node2tip.values())
    sd = statistics.stdev(node2tip.values())
    threshold = median + sd * maxdev
    to_keep = [k for k, v in node2tip.items() if v <= threshold]
    return to_keep


def analyse_leca(tr, leca, tname, leca_no, group2dom, max_sists=2):
    sisters_summary = analyse_sisters(tr, leca, group2dom, 2)
    rsl = leca.dist
    ebl = statistics.median([leca.get_distance(x) for x in leca.get_leaves()])
    if ebl == 0:
        ebl = statistics.mean([leca.get_distance(x) for x in leca.get_leaves()])
    
    if ebl == 0:
        return None

    basedict = {'Orthogroup': tname,
                'LECA': f'LECA_{leca_no}',
                'LECA_size': len(leca),
                'rsl': rsl,
                'ebl': ebl,
                'sl': rsl / ebl,
                'leca_is_root': leca.is_root(),
                'leca_support': leca.support,
                'leca_sister_support': leca.up.support}

    if sisters_summary is None:
        return None
    else:
        i = 0
        taxdict = {}
        for sister in sisters_summary:
            if i < 2:
                taxdict = {**basedict, **taxdict}
                for k, v in sisters_summary[sister].items():
                    if isinstance(v, dict):
                        dstr = ''
                        for m, n in v.items():
                            dstr += f'{m}: {n}\n'
                        taxdict[f'{sister}_{k}'] = dstr
                    else:
                        taxdict[f'{sister}_{k}'] = v

            i += 1

        leca_euks = [x.name for x in leca.get_leaves() if x.d == 'Eukaryota']
        criteria = is_LECA(leca_euks)
        clade = sisters_summary['S1']['MAB_group']
        odict = {**basedict,
                'S1_size': sisters_summary['S1']['size'],
                'donor': sisters_summary['S1']['mrca'],
                'S1_MAB_group_domain': sisters_summary['S1']['MAB_group_domain'],
                'S1_in_S2': sisters_summary['S1']['S1_in_S2'],
                **criteria,
                'LECA_seqs': ';'.join(leca.get_leaf_names()),
                'S1_seqs': sisters_summary['S1']['seqs']}

        return {'donor_information': odict, 'taxdict': taxdict}


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-t', '--trees', dest='trees',
                        help='File with multiple trees.')
    parser.add_argument('-l', '--lineage', dest='lineage',
                        help='Lineage file.')
    parser.add_argument('-p', '--prefix', dest='prefix',
                        help='Output prefix')
    args = parser.parse_args()

    trees = args.trees
    lng = args.lineage
    prefix = args.prefix

    create_folder(prefix.rsplit('/', 1)[0])
    lng = read_lineage(lng)
    group2dom = {v['sg']: v['d'] for k, v in lng.items()}

    odfl = []
    taxodfl = []
    otrees = ''
    # olecanodes = ''
    for tline in open(trees, 'r'):
        tline = read_treeline(tline)
        print(tline['tname'])
        tr = read_tree(tline['tree'], lng, get_sp=get_sp, mrca_annotation=False)

        doms = [x.d for x in tr.get_leaves()]
        if len(set(doms)) == 1:
            continue

        leca = get_leca(tr)

        i = 1
        while leca is not None:
            to_keep = [x for x in tr.get_leaf_names() if x not in leca[1].get_leaf_names()]
            if leca[0]:
                leca_analysis = analyse_leca(tr, leca[1], tline['tname'], i, group2dom, 2)
                if leca_analysis is not None:
                    odfl.append(leca_analysis['donor_information'])
                    taxodfl.append(leca_analysis['taxdict'])
                    otrees += f'{tline["tname"]}\tLECA_{i}\t{tline["model"]}\t{tr.write(features=["name"])}\n'
                    i += 1
            tr.prune(to_keep)
            leca = get_leca(tr)

    pd.DataFrame(odfl).to_csv(f'{prefix}_stats.tsv', sep='\t', index=False)
    pd.DataFrame(taxodfl).to_csv(f'{prefix}_taxonomy_S1.tsv', sep='\t', index=False)
    with open(f'{prefix}_LECA_trees.nwk', 'w') as f:
        f.write(otrees)

    return 0


if __name__ == '__main__':
    main()
