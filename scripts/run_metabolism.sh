for i in ../outputs/metabolism/*proteome.tsv
do
    bnm=$(basename $i "_KOori_proteome.tsv"); anvi-estimate-metabolism --enzymes-txt $i -O ../outputs/metabolism/${bnm} --include-kos-not-in-kofam
done
