#!/usr/bin/env python3

from glob import glob
import argparse as ap
import sys
from statistics import mean
import polars as pl

from utils import basename, create_folder
from sequence_utils import read_fasta, write_fasta, clean_seq


def trim_ends(sequence: str, non_masked_seq=None, soft=False):
    # Trim head
    seq = list(sequence)
    if seq[0] == 'X' or seq[1] == 'X':
        for i in range(0, len(seq)):
            if i != 0 and seq[i] != 'X':
                break
            start = i
        start += 1
    else:
        start = 0

    # Trim tail
    if seq[len(seq) - 1] == 'X' or seq[len(seq) - 2] == 'X':
        for i in range(len(seq) - 1, 0, -1):
            if i != len(seq) - 1 and seq[i] != 'X':
                break
            end = i
    else:
        end = len(seq)

    if non_masked_seq is None:
        return clean_seq(''.join(seq[start:end]))
    elif non_masked_seq is not None and soft:
        if len(non_masked_seq) != len(sequence):
            print('error')
        return clean_seq(''.join(non_masked_seq[start:end]))


def count_lowercase(string: str):
    return sum(1 for x in string if x.islower())


def has_stop_in(sequence: str):
    if '*' in sequence:
        seq = list(sequence)
        if seq.index('*') < len(seq) - 1:
            return True
        else:
            return False
    else:
        return False


def calculate_lcprop(sequence, soft=False):
    protlen = len(sequence)
    if soft:
        lc_count = count_lowercase(sequence)
    else:
        lc_count = list(sequence).count('X')

    if protlen > 0:
        lc_prop = lc_count / protlen
    else:
        lc_prop = None

    return {'lc_count': lc_count, 'lc_prop': lc_prop, 'protlen': protlen}


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-r', '--raw', dest='raw_proteomes',
                        help='Raw proteomes list.',
                        metavar='<path>',
                        nargs='*')
    parser.add_argument('-m', '--masked', dest='masked_proteomes',
                        help='Masked proteomes list with X in low complexity regions.',
                        metavar='<path>',
                        nargs='*')
    parser.add_argument('-o', '--output', dest='output',
                        help='Output file in tsv.',
                        metavar='<path.tsv>')
    parser.add_argument('-s', '--sequences', dest='sequences',
                        help='Output directory for the cleaned sequences.',
                        metavar='<path>')
    parser.add_argument('-f', '--soft', dest='soft',
                        help=('Write the cleaned proteome with the masked '
                              'amino acid, not X. By default it writes it with X.'),
                        action='store_true')
    parser.add_argument('-i', '--minimum', dest='minseqlen',
                        help=('Minimum sequence length'),
                        metavar='<N>',
                        type=int,
                        default=50)
    parser.add_argument('-a', '--maximum', dest='maxseqlen',
                        help=('Maximum sequence length'),
                        metavar='<N>',
                        type=int,
                        default=10000)
    parser.add_argument('-c', '--complexity', dest='complexity',
                        help=('Proportion of Xs in the proteome to consider it '
                              'a low-complexity protein'),
                        metavar='<N>',
                        type=float,
                        default=0.25)
    args = parser.parse_args()

    odir = args.sequences
    create_folder(odir)

    proteomes = {basename(x, '.fa'): x for x in args.raw_proteomes}
    masked = {basename(x, '.fa'): x for x in args.masked_proteomes}

    # Check that all the lists are the same
    if len(set(proteomes.keys()) & set(masked.keys())) != len(proteomes):
        sys.stdout('Not the same number of files or with different names.')
        sys.exit(1)

    # For each proteome
    odfl = []
    for proteome in proteomes:
        # Read raw and masked
        raw_seqs = read_fasta(proteomes[proteome], clean=False)
        masked_seqs = read_fasta(masked[proteome], clean=False)

        if not all([x in masked_seqs.keys() for x in raw_seqs.keys()]):
            sys.stdout('The proteins in the masked and raw genomes do not coincide.')
            sys.exit(1)

        oproteome = {}
        odfd = {'proteome': proteome,
                'proteome_length': len(raw_seqs),
                'raw_lc_proteins': 0,
                'masked_trimmed_lc_proteins': 0,
                'clean_proteome_length': 0}

        raw_lc_prop = []
        masked_trimmed_lc_prop = []

        raw_protlen = []
        masked_trimmed_protlen = []

        oproteome_protlen = []
        for protein in raw_seqs:
            if len(masked_seqs[protein]) == 1:
                continue
            # Trim ends of the masked protein
            print(masked_seqs[protein])
            trimmed_protein = trim_ends(masked_seqs[protein])
            trimmed_tofile = trim_ends(masked_seqs[protein], raw_seqs[protein],
                                       args.soft)

            # Calculating low complexity stats
            raw_lcstats = calculate_lcprop(raw_seqs[protein])
            masked_lcstats = calculate_lcprop(trimmed_protein)

            raw_lc_prop.append(raw_lcstats['lc_prop'])

            raw_protlen.append(raw_lcstats['protlen'])
            masked_trimmed_protlen.append(len(trimmed_protein))

            if (raw_lcstats['lc_prop'] is not None and
                raw_lcstats['lc_prop'] >= args.complexity):
                odfd['raw_lc_proteins'] += 1

            if masked_lcstats['lc_prop'] is not None:
                masked_trimmed_lc_prop.append(masked_lcstats['lc_prop'])

            if (masked_lcstats['lc_prop'] is not None and
                masked_lcstats['lc_prop'] >= args.complexity):
                odfd['masked_trimmed_lc_proteins'] += 1

            if (masked_lcstats['lc_prop'] is not None and
                masked_lcstats['lc_prop'] < args.complexity and
                masked_lcstats['protlen'] >= args.minseqlen and
                masked_lcstats['protlen'] <= args.maxseqlen and
                not has_stop_in(raw_seqs[protein])):
                oproteome[protein] = trimmed_tofile
                odfd['clean_proteome_length'] += 1
                oproteome_protlen.append(len(trimmed_protein))


        odfd['raw_lc_prots_prop'] = odfd['raw_lc_proteins'] / odfd['proteome_length']
        odfd['masked_lc_prots_prop'] = odfd['masked_trimmed_lc_proteins'] / odfd['proteome_length']
        odfd['lc_index'] = abs(odfd['masked_lc_prots_prop'] - odfd['raw_lc_prots_prop'])

        odfd['mean_raw_lc_prop'] = mean(raw_lc_prop)
        odfd['min_raw_lc_prop'] = min(raw_lc_prop)
        odfd['max_raw_lc_prop'] = max(raw_lc_prop)

        odfd['mean_masked_trimmed_lc_prop'] = mean(masked_trimmed_lc_prop)
        odfd['min_masked_trimmed_lc_prop'] = min(masked_trimmed_lc_prop)
        odfd['max_masked_trimmed_lc_prop'] = max(masked_trimmed_lc_prop)

        odfd['mean_raw_protlen'] = mean(raw_protlen)
        odfd['min_raw_protlen'] = min(raw_protlen)
        odfd['max_raw_protlen'] = max(raw_protlen)
        
        odfd['mean_masked_trimmed_protlen'] = mean(masked_trimmed_protlen)
        odfd['min_masked_trimmed_protlen'] = min(masked_trimmed_protlen)
        odfd['max_masked_trimmed_protlen'] = max(masked_trimmed_protlen)

        odfl.append(odfd)

        write_fasta(oproteome, f'{odir}/{proteome}.fa')
    
    odf = pl.DataFrame(odfl)
    odf.write_csv(args.output, separator='\t')

if __name__ == '__main__':
    main()