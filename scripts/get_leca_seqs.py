#!/usr/bin/env python3

from sequence_utils import read_fasta, write_fasta, get_seqs
from utils import basename, create_folder
from glob import glob

def get_sp(header):
    return header.split('_')[0]

create_folder(f'outputs/large_ogs/round_2/seqs/')

for db in ['TOLDBA', 'TOLDBB', 'TOLDBC']:
    gnms = {basename(x, '.fa'): x for x in glob('outputs/toldb/toldb_final_proteomes/*')}
    used_gnms = {}

    create_folder(f'outputs/large_ogs/round_2/seqs/{db}/')

    for line in open(f'outputs/large_ogs/round_1/{db}_r1_LECA_seqs.tsv', 'r'):
        og, seqs = line.strip().split('\t')
        ofile = og.replace('.final_LECA', '')
        odir = f'outputs/large_ogs/round_2/seqs/{db}/{ofile}'
        create_folder(odir)
        seqs = seqs.split(';')
        sp2gnm = {get_sp(x): read_fasta(gnms[get_sp(x)]) for x in seqs if 'Prok|' not in x}
        oseqs = {}
        for seq in seqs:
            if 'Prok|' not in seq:
                oseqs[seq] = sp2gnm[get_sp(seq)][seq]
        write_fasta(oseqs, outfile=f'{odir}/{ofile}.fa')
        print(f'{odir}/{ofile}.fa written')