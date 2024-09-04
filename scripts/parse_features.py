#!/usr/bin/env python3

import polars as pl

donors = []
for line in open('outputs/metabolism_vircleaned/donor_domain_KEGG_mapper_reconstruction.tsv'):
    leading_sp = len(line) - len(line.lstrip())
    line = line.strip()
    if leading_sp == 0:
        pathtype = line
    elif leading_sp == 4:
        pathsubtype = line
    elif leading_sp == 8 and line.startswith('0'):
        pathway = line.split(' ', 1)[1]
    elif leading_sp == 8 and line.startswith('K'):
        KO = line
    elif leading_sp == 12:
        donors_line = line.split(': ')
        if len(donors_line) > 1:
            donors.append({'type': pathtype,
                           'subtype': pathsubtype,
                           'pathway': pathway,
                           'KO': KO,
                           'donor': donors_line[0],
                           'donor_genes': donors_line[1]})

odf = pl.DataFrame(donors)
odf.write_csv(file='outputs/metabolism_vircleaned/donor_domain_KEGG_donors.tsv', separator='\t')