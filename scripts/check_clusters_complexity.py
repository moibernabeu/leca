#!/usr/bin/env python3

import argparse as ap
from glob import glob
import sys

from clean_proteomes import trim_ends, calculate_lcprop
from utils import basename
from sequence_utils import read_fasta, write_fasta
from cluster_utils import read_mmseqs_clusters

def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-s', '--sequences', dest='sequences',
                        help='Input sequences.',
                        metavar='<path>',
                        nargs='*')
    parser.add_argument('-m', '--masked', dest='masked',
                        help='Input masked sequences.',
                        metavar='<path>',
                        nargs='*')
    parser.add_argument('-c', '--clusters', dest='clusters',
                        help='Input clusters',
                        metavar='<path>',
                        nargs='*')
    args = parser.parse_args()

    args.sequences = glob('/gpfs/projects/bsc40/current/mgil/documents/leca_V3/outputs/cleaned_proteomes/proteomes/*')
    args.masked = glob('/gpfs/projects/bsc40/current/mgil/documents/leca_V3/outputs/masked_gnms/*')
    args.clusters = glob('/gpfs/projects/bsc40/current/mgil/documents/leca_V3/outputs/clustered_proteomes/proteomes/*/*_cluster.tsv')

    proteomes = {basename(x, '.fa'): x for x in args.sequences}
    masked = {basename(x, '.fa'): x for x in args.masked}
    cluster_files = {basename(x, '_cluster.tsv'): x for x in args.clusters}

    # Check that all the lists are the same
    if len(set(proteomes.keys()) & set(masked.keys())) != len(proteomes):
        sys.stdout('Not the same number of files or with different names.')
        sys.exit(1)

    for i, proteome in enumerate(proteomes.keys()):
        clean_seqs = read_fasta(proteomes[proteome], clean=False)
        masked_seqs = read_fasta(masked[proteome], clean=False)
        clusters = read_mmseqs_clusters(cluster_files[proteome])

        print(proteome, len(clean_seqs), len(masked_seqs))

        for cluster in clusters:
            miau = {'rep': cluster}
            for seq in clusters[cluster]:
                stats = calculate_lcprop(trim_ends(masked_seqs[seq]))
                miau[seq] = (stats['lc_prop'], stats['protlen'])
            print(miau)

    return 0


if __name__ == '__main__':
    main()