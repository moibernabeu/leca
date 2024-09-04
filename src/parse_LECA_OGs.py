#!/usr/bin/env python3

from cluster_utils import read_orthofinder_clusters, cluster_summary
from LECA_utils import is_LECA, get_sp
from tree_utils import read_lineage
from utils import create_folder

import polars as pl

import argparse as ap


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input',
                        help='Orthofinder Orthogroups.txt output file.',
                        metavar='<path>')
    parser.add_argument('-l', '--lineage', dest='lineage',
                        help='Lineage file.',
                        metavar='<path>')
    parser.add_argument('-o', '--odir', dest='odir',
                        help='Outputs directory.',
                        metavar='<path>')
    args = parser.parse_args()

    create_folder(args.odir)

    clusters = read_orthofinder_clusters(args.input)
    lng = read_lineage(args.lineage)

    odfl = []
    for clust in clusters:
        clustsumm = cluster_summary(clusters[clust], get_sp, og=clust,
                                    lng=lng, ranks=['sg', 'div'])
        criteria = is_LECA(clusters[clust], lng)
        
        odfl.append({**clustsumm, **criteria})

    df = pl.DataFrame(odfl)
    fdf = df.filter(pl.col('Two_stems') & pl.col('Five_sp') & pl.col('Three_supergroups'))

    df.write_csv(f'{args.odir}/orthogroups_summary.tsv', separator='\t')
    fdf.write_csv(f'{args.odir}/orthogroups_summary_filtered.tsv', separator='\t')

    return 0


if __name__ == '__main__':
    main()