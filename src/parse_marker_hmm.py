#!/usr/bin/env python3

import sys
import argparse as ap
import pandas as pd
from statistics import median
from utils import create_folder

from sequence_utils import read_fasta, write_fasta, get_seqs


def get_sp(header):
    return header.split('_')[0]


def read_hmmtab(hmmtabfile, evalthr=10, get_sp=get_sp, sp_max=5):
    # Opening the hmm file and iterating through lines
    search_dict = {}
    for line in open(hmmtabfile, 'r'):
        if '#' not in line:
            line = line.strip()
            line = ' '.join(line.split()).split(' ')
            if len(line) < 18:
                # Not enough fields
                # sys.exit(1)
                break
            
            # First hit dictionary
            if line[0] not in search_dict:
                hit = {}
                hit['target'] = line[0]
                hit['target_acc'] = line[1]
                hit['query'] = line[2]
                hit['query_acc'] = line[3]
                hit['evalue'] = float(line[4])
                hit['score'] = line[5]
                hit['bias'] = line[6]

                if hit['query'] not in search_dict:
                    search_dict[hit['query']] = {}
                    sp_count = {}

                if not get_sp(hit['target']) in sp_count.keys() and sp_max is not None:
                    sp_count[get_sp(hit['target'])] = 1
                else:
                    sp_count[get_sp(hit['target'])] += 1

                if sp_max is not None:
                    incorporate = sp_count[get_sp(hit['target'])] <= sp_max
                    print(incorporate, get_sp(hit['target']), sp_count[get_sp(hit['target'])], sp_max)
                else:
                    incorporate = True

                if hit['evalue'] <= evalthr and incorporate:
                    search_dict[hit['query']][hit['target']] = hit

    return search_dict


def hmm_summary(hmm_dict, splist=None, get_sp=get_sp):
    sp_in_hmm = [get_sp(x) for x in hmm_dict]

    if splist is None:
        splist = list(set(sp_in_hmm))

    spcov = len(set(sp_in_hmm)) / len(splist)

    spcount = {}
    for sp in splist:
        spcount[sp] = sp_in_hmm.count(sp)

    return {'summary': {'seqs_no': len(sp_in_hmm),
                        'sp_no': len(set(sp_in_hmm)),
                        'sp_list_no': len(splist),
                        'spcov': spcov,
                        'single_copy': len(sp_in_hmm) == len(set(sp_in_hmm)),
                        **spcount},
            'sp_in_hmm': list(set(sp_in_hmm))}


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', dest='hmmtab',
                        help='HMM tabular output.',
                        metavar='<path>')
    parser.add_argument('-s', '--sequences', dest='seqsfile',
                        help='Sequences of the genomes.',
                        metavar='<X>')
    parser.add_argument('-m', '--maxseqs', dest='maxseqs',
                        help='Maximum number of sequences to retrieve from each genome',
                        metavar='<N>',
                        type=int)
    parser.add_argument('-l', '--splist', dest='splist',
                        help='Species list',
                        metavar='<path>')
    parser.add_argument('-t', '--tabout', dest='tabout',
                        help='Output path for the tsv of the presence absence.',
                        metavar='<path>')
    parser.add_argument('-o', '--seqsdir', dest='seqsdir',
                        help='Output directory for the sequences.',
                        metavar='<path>')
    args = parser.parse_args()

    splist = [x.split('\t')[1] for x in open(args.splist, 'r')]
    create_folder(args.seqsdir)

    hmm = read_hmmtab(args.hmmtab, sp_max=args.maxseqs)
    seqs = read_fasta(args.seqsfile, clean=True)

    odfl = []
    for marker in hmm:
        summ = hmm_summary(hmm[marker], splist)
        # Filtering by the species coverage of the found homologs
        if summ['summary']['spcov'] >= 0.8:
            byevalue = [(k, v['evalue']) for k, v in hmm[marker].items()]
            byevalue.sort(key=lambda x: x[1], reverse=False)
            marker_seqs = get_seqs(seqs, hmm[marker].keys())
            seq_lens = [len(marker_seqs[x]) for x in marker_seqs]
            oseqs = {}
            ospp = set()
            for seq, evalue in byevalue:
                # Filtering by the length of the sequences
                if (len(marker_seqs[seq]) >= 0.5 * median(seq_lens) and
                    len(marker_seqs[seq]) <= 2 * median(seq_lens)):
                    oseqs[seq] = marker_seqs[seq]
                    ospp.add(get_sp(seq))
                else:
                    print(f'{seq} has been removed due to length.')
                
                if len(set(ospp)) == len(summ['sp_in_hmm']):
                    print('Reached the species number')
                    break

            oseqs_sp = [get_sp(x) for x in oseqs]
            summ['filtered_spcov'] = len(set(oseqs_sp)) / len(splist)
            summ['included'] = summ['filtered_spcov'] >= 0.8
            if summ['included']:
                ofile = f'{args.seqsdir}/{marker}.fa'
                write_fasta(oseqs, ofile)

        odfl.append(summ['summary'])

    pd.DataFrame(odfl).to_csv(args.tabout, sep='\t', index=False)

    return 0


if __name__ == '__main__':
    main()
