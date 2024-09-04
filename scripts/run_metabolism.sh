for i in ../outputs/metabolism_vircleaned/*.tsv
do
    bnm=$(basename $i "_KOori_proteome.tsv"); anvi-estimate-metabolism --enzymes-txt $i -O ../outputs/metabolism_vircleaned/${bnm} --include-kos-not-in-kofam
done
