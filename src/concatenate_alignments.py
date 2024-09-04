#!/usr/bin/env python3

from sequence_utils import read_alignment, write_fasta
from utils import create_folder

from glob import glob
import argparse as ap

def get_sp(seqid):
    return seqid.split('_')[0]


def complete_alignment(seqs, sp_list):
    sp_in_seqs = [get_sp(x) for x in seqs]
    aln_length = len(list(seqs.values())[0])
    oseqs = {}

    for sp in sp_list:
        if sp not in sp_in_seqs:
            oseqs[sp] = '-' * aln_length

    for seq in seqs:
        sp = get_sp(seq)
        if sp in sp_in_seqs:
            oseqs[sp] = seqs[seq]

    return oseqs


def cat_alignments(alnlist, get_sp=get_sp):
    alns = {}
    for file in alnlist:
        filename = file.rsplit('/', 1)[1].rsplit('.', 1)[0]
        alns[filename] = read_alignment(file, True, get_sp)

    species = list(set([get_sp(y) for y in sum([list(alns[x].keys()) for x in alns], [])]))

    cat_alns = {}
    for aln in alns:
        cat_alns[aln] = complete_alignment(alns[aln], species)

    oseqs = {}
    summary = {}
    initial_position = 1
    for aln in cat_alns:
        for sp in species:
            if sp not in oseqs.keys():
                oseqs[sp] = cat_alns[aln][sp]
            else:
                oseqs[sp] += cat_alns[aln][sp]

        final_position = len(oseqs[sp])

        summary[aln] = (initial_position, final_position)

        print(aln, initial_position, final_position, sep='\t')

        initial_position = final_position + 1

    return {'cat_aln': oseqs, 'summary': summary}


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input',
                        help='Input alignments listed or as a wildcard.',
                        nargs='*')
    parser.add_argument('-o', '--output', dest='output',
                        help='Output concatenated alignment.')
    args = parser.parse_args()

    input_alns = args.input
    output = args.output

    create_folder(output.rsplit('/', 1)[0])

    cat_aln = cat_alignments(alnlist=input_alns, get_sp=get_sp)
    write_fasta(cat_aln['cat_aln'], output)

    return 0

if __name__ == '__main__':
    main()