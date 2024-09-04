#!/usr/bin/env python3

from tree_utils import read_lineage, read_tree
from LECA_utils import get_sp
from utils import get_list_props, count_list_items

import ete3
import sys
from statistics import mean
import polars
import argparse as ap

def family_name(x):
    return x.split('/')[13]

def read_leca_treeline(line):
    family, leca, model, nwk = line.strip().split('\t')
    odict = {'family': family, 'LECA': leca, 'model': model, 'nwk': nwk}
    return odict

def get_outgroup(tr: ete3.PhyloTree):
    og_seqs = tr.get_tree_root().children[0].get_leaf_names()
    return og_seqs

def set_ml_outgroup(tr: ete3.PhyloTree, og_seqs):
    mrca = tr.get_common_ancestor(og_seqs)
    if len(mrca) > 2 * len(og_seqs):
        tr.set_outgroup(og_seqs[0])
    else:
        tr.set_outgroup(mrca)
    return tr

def get_leca_split(tr: ete3.PhyloTree):
    lecaup = tr.search_nodes(name='LECA')[0].up
    if lecaup is not None:
        split = lecaup.children
        split_dict = {}
        for x in split:
            if x.name == '' or len(x) == 1:
                split_dict['S1'] = x
            elif x.name == 'LECA':
                split_dict['LECA'] = x
        return split_dict
    else:
        return None

def annotate_btr_leca(btr: ete3.PhyloTree, leca_seqs):
    btr.get_common_ancestor(leca_seqs).add_feature('name', 'LECA')
    return 0

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def jaccard(lst1, lst2):
    inters = intersection(lst1, lst2)
    uni = list(set(lst1 + lst2))
    return len(inters) / len(uni)

def analyse_sister(sister, group2dom):
    domains = [x.d for x in sister.get_leaves()]

    domains = get_list_props(domains, sort=True)

    groups = [x.sg for x in sister.get_leaves() if x.d != 'Eukaryota']
    if len(groups) > 0:
        groups = get_list_props(groups, sort=True)

        try:
            mixed = groups[0][1] == groups[1][1]
        except IndexError:
            mixed = False

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
                 'support': sister.support}

        return sdict
    else:
        return None


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-t', '--trees', dest='mltrees',
                        help='ML trees')
    parser.add_argument('-b', '--bstrees', dest='bstrees',
                        help='Bootstrap trees files list.')
    parser.add_argument('-o', '--output', dest='ofile',
                        help='Output table, in tsv.')
    args = parser.parse_args()

    lngfile = '/gpfs/projects/bsc40/current/mgil/documents/leca_V3/data/LECA_proteophylum.lng'
    mltrees = args.mltrees
    bootstrap_files = args.bstrees
    ofile = args.ofile

    # mltrees = '../outputs/final_trees/TOLDBA_LECA_trees.nwk'
    # bootstrap_files = '../outputs/taxonomic_bootstrap/TOLDBA_ufboot_files.txt'
    # ofile = '../outputs/taxonomic_bootstrap/TOLDBA.tsv'
    
    bootstrap_files = {family_name(x.strip()): x.strip() for x in open(bootstrap_files, 'r')}
    lng = read_lineage(lngfile)
    group2dom = {v['sg']: v['d'] for k, v in lng.items()}

    odfl = []
    for i, line in enumerate(open(mltrees, 'r')):
        treeline = read_leca_treeline(line)
        print(f'Parsing tree {i}: {treeline["family"]}, {treeline["LECA"]}')
        tr = read_tree(treeline['nwk'], lng, get_sp=get_sp, root=False, mrca_annotation=False)
        og = get_outgroup(tr)
        ml_split = get_leca_split(tr)
        leca_seqs = [x.name for x in ml_split['LECA'].get_leaves() if x.d == 'Eukaryota']

        ml_sister = analyse_sister(ml_split['S1'], group2dom)

        boottrees = {'total': 0, 'S1_jaccard': [], 'LECA_jaccard': [],
                    'donor': [], 'donor_domain': [], 'LECA_seqs': [], 'LECA_bs': 0}
        for boottree in open(bootstrap_files[treeline['family']], 'r'):
            boottree = boottree.strip()
            btr = read_tree(boottree, lng, get_sp=get_sp, root=False, mrca_annotation=False)
            set_ml_outgroup(btr, og)
            annotate_btr_leca(btr, leca_seqs)
            btr_split = get_leca_split(btr)
            if btr_split is not None:
                sister_analysis = analyse_sister(btr_split['S1'], group2dom)
                LECA_jaccard = jaccard(ml_split['LECA'].get_leaf_names(),
                                        btr_split['LECA'].get_leaf_names())
                if sister_analysis is not None:
                    S1_jaccard = jaccard(ml_split['S1'].get_leaf_names(),
                                        btr_split['S1'].get_leaf_names())
                    boottrees['S1_jaccard'].append(S1_jaccard)
                    boottrees['LECA_jaccard'].append(LECA_jaccard)
                    boottrees['donor'].append(sister_analysis['mrca'])
                    boottrees['donor_domain'].append(sister_analysis['MAB_group_domain'])
                    boottrees['LECA_seqs'].append(btr_split['LECA'].get_leaf_names())
                    boottrees['total'] += 1
            if LECA_jaccard == 1:
                boottrees['LECA_bs'] += 1

        boottrees['LECA_bs'] = boottrees['LECA_bs'] / 1000

        donor_props = get_list_props(boottrees['donor'])
        if ml_sister['mrca'] in donor_props.keys():
            ml_donor_bs = donor_props[ml_sister['mrca']]
        else:
            ml_donor_bs = 0
    
        if 'Asgardarchaeota' in donor_props:
            asgard = donor_props['Asgardarchaeota']
        else:
            asgard = 'NA'
        if 'Thermoproteota' in donor_props:
            thermo = donor_props['Thermoproteota']
        else:
            thermo = 'NA'
        if 'Alphaproteobacteria' in donor_props:
            alpha = donor_props['Alphaproteobacteria']
        else:
            alpha = 'NA'
        if 'Gammaproteobacteria' in donor_props:
            gamma = donor_props['Gammaproteobacteria']
        else:
            gamma = 'NA'
        if 'Myxococcota' in donor_props:
            myxo = donor_props['Myxococcota']
        else:
            myxo = 'NA'
        if 'Nucleocytoviricota' in donor_props:
            nucl = donor_props['Nucleocytoviricota']
        else:
            nucl = 'NA'
        if 'Planctomycetota' in donor_props:
            planct = donor_props['Planctomycetota']
        else:
            planct = 'NA'
        if 'Acidobacteriota' in donor_props:
            acido = donor_props['Acidobacteriota']
        else:
            acido = 'NA'
        if 'Chloroflexota' in donor_props:
            chloro = donor_props['Chloroflexota']
        else:
            chloro = 'NA'

        donor_props = sorted(donor_props.items(), key=lambda kv: kv[1], reverse=True)

        btr_LECA_seqs = count_list_items(sum(boottrees['LECA_seqs'], []))
        btr_LECA_core = [k for k, v in btr_LECA_seqs.items() if v / 1000 >= 0.25]

        odfd = {'family': treeline['family'],
                'LECA': treeline['LECA'],
                'LECA_bs': ml_split['LECA'].support,
                'LECA_bs_us': boottrees['LECA_bs'],
                'used_btr': boottrees['total'] / 1000,
                'S1_jaccard_mean': mean(boottrees['S1_jaccard']),
                'LECA_jaccard_mean': mean(boottrees['LECA_jaccard']),
                'ML_donor': ml_sister['mrca'],
                'btr_ML_donor_prop': ml_donor_bs,
                'btr_majoritary_donor': donor_props[0][0],
                'btr_majoritary_donor_prop': donor_props[0][1],
                'asgard': asgard,
                'alpha': alpha,
                'gamma': gamma,
                'thermo': thermo,
                'myxo': myxo,
                'nucl': nucl,
                'planct': planct,
                'acido': acido,
                'chloro': chloro,
                'btr_LECA_seqs_no': len(btr_LECA_seqs),
                'btr_LECA_core_no': len(btr_LECA_core),
                'LECA_core': ';'.join(btr_LECA_core)}

        odfl.append(odfd)

    df = polars.DataFrame(odfl, strict=False)
    df.write_csv(ofile, separator='\t')

    return 0


if __name__ == '__main__':
    main()
