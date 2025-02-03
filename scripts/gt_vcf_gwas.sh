#!/bin/bash

# Assign command-line arguments to variables
inds=$1.inds
chr=$2
midchrom=$3
gen1=$4
gen2=$5

## if ${1} contains "_gwas" then proceed
# This section for GWAS
# Creates phenotype file for samples
## Else if iHS then just remove .pos, .inds, .gt

if [[ "$1" == *"gwas"* ]]; then
	# mv ${1}.vcf.gz ${1}_${4}_${5}.vcf.gz
	# bcftools index -f ${1}_${4}_${5}.vcf.gz
	 sort -k1n ${1}.inds | mawk '{ if ($2 == "Dead") {print $1"\t"$1"\t"$2"\t1"}  else {print $1"\t"$1"\t"$2"\t2"} }' > ${1}_${4}_${5}.pheno
	 #rm ${1}.gt
	 rm ${1}.inds

    plink --vcf ${1}_${4}_${5}.vcf.gz \
    --allow-no-sex \
    --geno 0.001 \
    --maf 0.05 \
    --logistic \
    --adjust \
    --seed 62 \
    --pheno ${1}_${4}_${5}.pheno \
    --mpheno 2 \
    --set-missing-var-ids @_# \
    --out ${1}_add

    grep -v NA ${1}_add.assoc.logistic | mawk 'NR==FNR{c[$2]++;next};c[$2]' - ${1}_add.assoc.logistic.adjusted | sort -n -t_ -k1,1 -k2,2 > ${1}_add.tmp
    grep -v NA ${1}_add.assoc.logistic | paste - ${1}_add.tmp | tr -s ' ' '\t' | cut -f1,11,12,13 --complement | gzip > ${1}_${4}_${5}_add.assoc.out.gz
    rm ${1}_add.tmp
    rm ${1}_add.assoc.logistic
    rm ${1}_add.assoc.logistic.adjusted


    plink --vcf ${1}_${4}_${5}.vcf.gz \
    --allow-no-sex \
    --geno 0.001 \
    --maf 0.05 \
    --logistic recessive \
    --model rec \
    --adjust \
    --seed 62 \
    --pheno ${1}_${4}_${5}.pheno \
    --mpheno 2 \
    --set-missing-var-ids @_# \
    --out ${1}_rec

    rm *.nosex

    grep -v NA ${1}_rec.assoc.logistic | mawk 'NR==FNR{c[$2]++;next};c[$2]' - ${1}_rec.assoc.logistic.adjusted | sort -n -t_ -k1,1 -k2,2 > ${1}_rec.tmp
    grep -v NA ${1}_rec.assoc.logistic | paste - ${1}_rec.tmp | tr -s ' ' '\t' | cut -f1,11,12,13 --complement | gzip > ${1}_${4}_${5}_rec.assoc.out.gz
    rm ${1}_rec.tmp
    rm ${1}_rec.assoc.logistic
    rm ${1}_rec.assoc.logistic.adjusted
    gzip --force ${1}_rec.model
fi
wait
