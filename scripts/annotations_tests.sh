for i in outputs/mLECA_origins/TOLDB*_stats.tsv;do bnm=$(basename $i "_stats.tsv"); awk -F '\t' '{if (NR != 1) {print $1 "_" $2 "\t" $24}}' $i > outputs/functional_annotation/${bnm}_mLECAOGs_sister.tsv;done

for i in outputs/mLECA_origins/TOLDB*_stats.tsv;do bnm=$(basename $i "_stats.tsv"); awk -F '\t' '{if (NR != 1) {print $1 "_" $2 "\t" $23}}' $i > outputs/functional_annotation/${bnm}_mLECAOGs.tsv;done

for i in outputs/functional_annotation/TOLDB*_mLECAOGs.tsv;do bnm=$(basename $i ".tsv"); python3 src/annotate_LECA_OGs.py -a /gpfs/projects/bsc40/current/mmarcet/LECA3/metabolism/predictions_KO_unique_euk.txt -g $i -o outputs/functional_annotation/${bnm}_annotation.tsv -e;done

for i in outputs/functional_annotation/TOLDB*_mLECAOGs_sister.tsv;do bnm=$(basename $i ".tsv"); python3 src/annotate_LECA_OGs.py -a /gpfs/projects/bsc40/current/mgil/documents/leca_V3/outputs/functional_annotation/prok_KOs.tsv -g $i -o outputs/functional_annotation/${bnm}_annotation.tsv;done

for i in outputs/functional_annotation/TOLDB*_mLECAOGs.tsv;do bnm=$(basename $i ".tsv"); python3 src/annotate_LECA_OGs.py -a /gpfs/projects/bsc40/current/mgil/documents/leca_V3/outputs/cog_annotation/proteins2COG.tsv -g $i -o outputs/functional_annotation/${bnm}_annotation_COG.tsv -e;done

# To run in tests/innovations_discussion
cut -f 1 /gpfs/projects/bsc40/shared_projects/LECA/V3/02.LECA_OGs/TOLDBC.LECAgroups.txt > 01_TOLDBC_LECA_OGs.txt

for i in ../../../../../../shared_projects/LECA/V3/initial_trees/trees/TOLDBC/Folder_*/*;do basename $i;done > 02_TOLDBC_LECA_OGs_to_expand.txt

for i in ../../../../../../shared_projects/LECA/V3/initial_trees/trees/TOLDBC/Folder_*/*/*.final.fa;do basename $i ".final.fa";done > 03_TOLDBC_LECA_OGs_expanded.txt

for i in ../../../../../../shared_projects/LECA/V3/initial_trees/trees/TOLDBC/Folder_*/*/*/*/*final.nwk;do basename $i ".final.fa";done > 04_TOLDBC_LECA_OGs_expanded_trees.txt

cut -f 1 ../../outputs/mLECA_origins/TOLDBC_stats.tsv | cut -f 1 -d '_' | awk 'NR != 1' | sort -u > 05_TOLDBC_acquisitions.txt

for i in outputs/metabolism/TOLDBA_3sg_KOori_proteome.tsv -O outputs/metabolism/*;do bnm=$(basename $i "_KOori_proteome.tsv"); anvi-estimate-metabolism --enzymes-txt $i -O outputs/metabolism/${bnm} --include-kos-not-in-kofam; done
