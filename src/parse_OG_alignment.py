#!/usr/bin/env python3

from sequence_utils import read_alignment, write_fasta, get_seqs
from LECA_utils import is_LECA
from utils import basename, create_folder

from statistics import quantiles
import argparse as ap
import polars as pl

canonic_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
              'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def seq_completeness(sequence: str, unambiguous: list=canonic_aa):
    sequence = list(sequence)

    ambiguous = 0
    for aa in sequence:
        if aa not in unambiguous:
            ambiguous += 1

    return 1 - (ambiguous / len(sequence))

def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input',
                        help='Alignment.',
                        metavar='<path>',
                        nargs='*')
    parser.add_argument('-o', '--output', dest='ofile',
                        help='Output summary table.',
                        metavar='<path.tsv>')
    parser.add_argument('-s', '--seqsdir', dest='seqsdir',
                        help='Output directory for the clean sequences',
                        metavar='<path>')
    args = parser.parse_args()
 
    create_folder(args.seqsdir)

    odfl = []
    for alnfile in args.input:
        print(f'Parsing {alnfile}')
        og = basename(alnfile, '.aln')
        aln = read_alignment(alnfile)

        seqs_compl = {}
        for seq in aln:
            seqs_compl[seq] = seq_completeness(aln[seq])

        to_keep = []
        compl_thr = quantiles(seqs_compl.values(), n = 10)[0]
        for seq in seqs_compl:
            if seqs_compl[seq] >= compl_thr or seqs_compl[seq] >= 0.75:
                to_keep.append(seq)

        criteria = is_LECA(to_keep)
        odfd = {'og': og,
                'file': alnfile,
                'pre_aln_seqsno': len(aln),
                **{f'pre_{k}': v for k, v in is_LECA(aln.keys()).items()},
                'post_aln_seqsno': len(to_keep),
                **{f'post_{k}': v for k, v in criteria.items()},
                'to_keep': ';'.join(to_keep)}

        if criteria['is_LECA'] and args.seqsdir is not None:
            oseqs = {k: v.replace('-', '') for k, v in aln.items()}
            write_fasta(get_seqs(oseqs, to_keep), f'{args.seqsdir}/{og}.fa')

        odfl.append(odfd)

    df = pl.DataFrame(odfl)
    df.write_csv(args.ofile, separator='\t')

    return 0


if __name__ == '__main__':
    main()
