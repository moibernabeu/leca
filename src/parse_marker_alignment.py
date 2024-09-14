#!/usr/bin/env python3

#!/usr/bin/env python3

from sequence_utils import read_alignment, write_fasta, get_seqs
from LECA_utils import is_LECA, get_sp
from utils import basename, create_folder
from pprint import pprint

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


def solitary_seqs(aln: dict, unambiguous: list=canonic_aa,
                  completeness_threshold: float=0.75):
    # Getting the alignment length
    aln_len = len(aln[list(aln.keys())[0]])

    # Retrieving the unambiguous positions
    unambiguous_positions = {}
    for i in range(0, aln_len):
        position = list()
        seqs = []
        for seq in aln:
            position.append(aln[seq][i])
            if aln[seq][i] not in unambiguous:
                seqs.append(seq)

        # Checking whether the position is ambiguous
        if 1 - seq_completeness(position) >= completeness_threshold:
            # Counting the non ambiguous sequences in the ambiguous site
            # (solitary site)
            for seq in seqs:
                if seq not in unambiguous_positions:
                    unambiguous_positions[seq] = 0
                unambiguous_positions[seq] += 1

    # Calculating the non-solitary sites proportion per each sequences
    solitary_positions_prop = {}
    for seq in unambiguous_positions:
        solitary_positions_prop[seq] = unambiguous_positions[seq] / aln_len

    solitary_threshold = quantiles(solitary_positions_prop.values(), n = 8)[0]

    # Returning the list of sequences with a percentage on non-solitary
    # positions lower than the calculated threshold. They are the sequence that
    # are generating gaps, as they have many non-ambiguous characters
    # in ambiguous positions.
    return [k for k, v in solitary_positions_prop.items() if v <= solitary_threshold]

def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input',
                        help='Input alignments.',
                        nargs='*')
    parser.add_argument('-s', '--seqsdir', dest='oseqs',
                        help='Output directory for the sequences.')
    parser.add_argument('-o', '--output', dest='output',
                        help='Output table in .tsv format.')
    parser.add_argument('-l', '--splist', dest='splist',
                        help='Species list.')
    args = parser.parse_args()

    create_folder(args.oseqs)
    splist = [x.split('\t')[1] for x in open(args.splist, 'r')]

    odfl = []
    for alignment in args.input:
        aln = read_alignment(alignment, get_sp=get_sp)

        # Calculating sequence completeness
        seqs_compl = {}
        for seq in aln:
            seqs_compl[seq] = seq_completeness(aln[seq])

        to_keep = []
        compl_thr = quantiles(seqs_compl.values(), n = 10)[0]
        for seq in seqs_compl:
            # Adding the sequences that have a completeness higher than the
            # calculated threshold or higher than 0.75 to the list of sequences
            # to keep
            if seqs_compl[seq] >= compl_thr or seqs_compl[seq] >= 0.75:
                to_keep.append(seq)

        # Analysing the sequences to discard as they add gaps to the alignment.
        to_discard = solitary_seqs(aln, completeness_threshold=compl_thr)
        to_keep = [x for x in to_keep if x not in to_discard]

        # Writing the output filtered sequences
        oseqs = {}
        for seq in aln:
            if seq in to_keep:
                oseqs[seq] = aln[seq].replace('-', '')

        aln_spp = [get_sp(x) for x in aln]
        out_spp = [get_sp(x) for x in aln]
        spp_count = {}
        for spp in splist:
            spp_count[spp] = out_spp.count(spp)
            
        bnm = basename(alignment, '.aln')
        ofile = f'{args.oseqs}/{bnm}.fa'

        odict = {}
        odict['gene'] = bnm
        odict['aln_len'] = len(aln)
        odict['cleaned_len'] = len(oseqs)
        odict['removed_seqs'] = len(aln) - len(oseqs)
        odict['aln_spp_no'] = len(set(aln_spp))
        odict['cleaned_spp_no'] = len(set(out_spp))
        odict['aln_spp_prop'] = len(set(aln_spp)) / 100
        odict['aln_spp_no'] = len(set(out_spp)) / 100
        odict = {**odict, **spp_count}
        odict['cleaned_seqs_file'] = ofile
        odfl.append(odict)

        write_fasta(oseqs, outfile=ofile)
        print(f'{ofile} written')

    odf = pl.from_dicts(odfl)
    odf.write_csv(file=args.output, separator='\t')

    return 0


if __name__ == '__main__':
    main()