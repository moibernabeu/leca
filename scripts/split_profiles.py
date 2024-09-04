#!/usr/bin/env python3

odir = '/gpfs/projects/bsc40/databases/COG_profiles'

currprof = ''
i = 0
for line in open('/gpfs/projects/bsc40/current/mgil/dbs/profiles/NCBI_cogs.hmm', 'r'):
    if line.startswith('NAME'):
        cogname = " ".join(line.strip().split()).split(' ', 1)[1]

    if line.strip() == '//':
        currprof += line
        print(f'Writing COG {cogname}')
        with open(f'{odir}/{cogname}.hmm', 'w') as f:
            f.write(currprof)
        i += 1
        currprof = ''
    else:
        currprof += line

    # if i == 10:
    #     break