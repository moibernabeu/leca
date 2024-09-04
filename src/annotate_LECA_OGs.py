#!/usr/bin/env python3

'''
The format of the KOs and orthogroups must be:
protein_header<tab>KO1;KO2;...
orthogroup<tab>protein1;protein2;...
'''

import polars as pl
import argparse as ap

from utils import create_folder


def rank_score(within_props: list, background_props: dict):
    # Calculate ratio
    ratio = []
    for id, prop in within_props:
        ratio.append((id, prop / background_props[id]))

    # Convert the ratio to a dataframe, sorting and obtaining the index
    ratiodf = pl.DataFrame(ratio, schema={'id': pl.String,
                                          'ratio': pl.Float32})
    ratiodf = ratiodf.sort(by=['ratio'], descending=True)
    ratiodf = ratiodf.with_row_index(name='rank_ratio')

    # Convert the proportions to a dataframe, sorting and obtaining the index
    propsdf = pl.DataFrame(within_props, schema={'id': pl.String,
                                                 'prop': pl.Float32})
    propsdf = propsdf.sort(by=['prop'], descending=True)
    propsdf = propsdf.with_row_index(name='rank_prop')

    # Join the ratio and proportions dataframe
    joined = ratiodf.join(propsdf, on='id')

    # Obtaining the score, that is the sum of the ranks in
    # the proportion and the ratio
    joined = joined.with_columns(score=pl.col('rank_prop') + pl.col('rank_ratio'))
    joined = joined.sort(by='score', descending=False)

    if len(joined) > 0:
        # Filtering the joined dataframe and 
        assigned = list(joined.filter(pl.col('score') == min(joined['score']))['id'])

        return assigned
    else:
        return None


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-a', '--annotations', dest='annotations',
                        help=('Annotations file, output from the merger of '
                              'KOfamscan and eggNOG results.'))
    parser.add_argument('-g', '--orthogroups', dest='orthogroups',
                        help='Orthofinder orthogroups parser output.')
    parser.add_argument('-o', '--output', dest='output',
                        help='Output file in tsv.')
    parser.add_argument('-e', '--onlyeuk', dest='onlyeuk',
                        help='Whether to remove prokaryotic sequences form the orthogroup.',
                        action='store_true')
    args = parser.parse_args()

    create_folder(args.output.rsplit('/', 1)[0])

    prots_KOs = {}
    all_kos = []
    background_count = {}
    for line in open(args.annotations, 'r'):
        protein, KOs = line.strip().split('\t', 1)
        KOs = KOs.split(';')

        if len(KOs) > 0:
            prots_KOs[protein] = KOs
            all_kos += KOs
            for ko in KOs:
                if ko not in background_count.keys():
                    background_count[ko] = 0
                background_count[ko] += 1

    print('Prots KOs retrieved')

    total_kos = sum(background_count.values())

    background_props = {k: v / total_kos for k, v in background_count.items()}

    print('Props calculated')

    dfl = []
    # Iterating through orthogroups
    for line in open(args.orthogroups, 'r'):
        og, seqs = line.strip().split('\t', 1)
        seqs = seqs.split(';')
        if args.onlyeuk:
            seqs = [x for x in seqs if 'Prok' not in x]
        dfd = {'og': og}
        print(f'Analysing OG: {og}')
        # Extracting the KOs from the OG
        within_KOs = {}
        prots_without_KOs = 0
        for prot in seqs:
            if prot in prots_KOs.keys():
                for ko in prots_KOs[prot]:
                    if ko not in within_KOs.keys():
                        within_KOs[ko] = 0
                    within_KOs[ko] += 1
            else:
                prots_without_KOs += 1

        print('KOs extracted, calculating proportions')

        # Calculating the KO proportions within the groups
        total_within_kos = sum(within_KOs.values())
        within_props = [(k, v / total_within_kos) for k, v in within_KOs.items()]

        print('Calculating the assigned KO.')

        # Getting the assigned KOs
        assigned_kos = rank_score(within_props, background_props)

        if assigned_kos is not None:
            KOs_str = ';'.join(assigned_kos)
            print(assigned_kos, [x[1] for x in within_props if x[0] in assigned_kos])
            KOs_prop = [x[1] for x in within_props if x[0] in assigned_kos]
            prop = max(KOs_prop)
            dfd = {'Orthogroup': og,
                   'KOs': KOs_str,
                   'prots_no': len(seqs),
                   'prots_with_KOs': float(len(seqs) - prots_without_KOs),
                   'pop_annotated': (len(seqs) - prots_without_KOs) / len(seqs),
                   'Max_freq_assigned_KO': float(prop)}
        else:
            dfd = {'Orthogroup': og,
                   'KOs': 'NA',
                   'prots_no': len(seqs),
                   'prots_with_KOs': float('nan'),
                   'pop_annotated': (len(seqs) - prots_without_KOs) / len(seqs),
                   'Max_freq_assigned_KO': float('nan')}

        dfl.append(dfd)

    df = pl.DataFrame(dfl)
    df.write_csv(args.output, separator='\t')

    return 0


if __name__ == '__main__':
    main()
