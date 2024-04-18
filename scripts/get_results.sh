#!/bin/bash

echo -e "\nepiSelect: Get results for multiple simulations and resamplings.\n"
# usage-message
: ${1?"Usage: $0 <model> <generation1> <generation2> <sample_size> [demeGeneration] [midchrom] [nsims] [nresamp] [workdir]"
"NOTE: If different sample sizes per generation/group, provide '_' between values "
"e.g. bash $0 add_grad_v1.0_f0.1_sim 40001 40003 10_15"
}

if [[ $# -lt 4 && $# -gt 0 ]] ; then
    echo -e 'You are missing arguments.\nUsage: bash power_calculate.sh <model> <generation1> <generation2> <sample_size> [demeGeneration] [midchrom] [nsims] [nresamp] [workdir]'
    exit 0
elif [[ $# -gt 4 && $# -lt 9 ]] ; then
    echo 'You can either provide all 9 arguments OR the first 4 with adjusted parameters within the script.'
    exit 0
elif [[ $# -eq 4 ]] ; then
    ## CHANGE THESE PARAMETERS FOR USER NEEDS
    deme=40000
    midchrom=33759515
    nsims=10
    nresamp=100
    workd="/home/lnd775/data/epiSelect/plague/results/"
else
    deme=$5
    midchrom=$6
    nsims=$7
    nresamp=$8
    workd=$9        
fi

model=$1
gen1=$2
gen2=$3
# If more than one sample size, provide "_" between values
sample=$4   
# Get inheritance mode
ih=$(echo ${model} | cut -d'_' -f1)
# Get viability and adjust for directory name
VAA=$(echo ${model} | cut -d'_' -f3 | cut -d'v' -f2)
if [[ ${VAA} == "1.0" ]] ; then
    VAA=1
else
    VAA=$(echo ${VAA} | tr -d ".")
fi

# Where we output the power calculations
plots=${workd}"power/"${ih}"/VAA"${VAA}


if [[ ${sample} == *"_"* ]]; then
	samples=${sample}
	sample1=$(cut -d'_' -f1 <<< ${sample})
	sample2=$(cut -d'_' -f2 <<< ${sample})
	sample=$(echo ${sample1}"_n"${sample2})
else
	samples=${sample}
	sample=$(echo ${sample}"_n"${sample})
fi

cd $plots/
if [ ! -d "${gen1}_${gen2}" ]; then
    mkdir ${gen1}_${gen2}
fi

for ((i=1;i<=nsims;i++)); do
sim=${i}
gravel=$sim
for ((k=1;k<=resamp; k++)); do
selected=$(sed -n 1p ${workd}${model}${i}/selected.txt)
freq=$(grep -w ${deme} ${workd}${model}${i}/${model}${i}.out | sed -n 2p | awk '{print $3}')

if [[ ${selected} -ge ${midchrom} ]]; then
    selected=$(echo ${selected}-${midchrom} | bc -l)
    selected="chr22_"${selected}
    chr=22
else
    selected="chr21_"${selected}
    chr=21
fi

# FST results
tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/fst_${gen1}_${gen2}_top_cand.txt | awk -v k="$k" -v gravel="$gravel" -v selected="$selected" -v freq="$freq" '{print $1"\t"$2"\t"$7"\t"$3"\t"NR"\t"gravel"\t"selected"\tFST\tresamp_"k"\t"freq"\t"$6}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topTen_mastertable.txt
# JSFS results
tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/jsfs_${gen1}_${gen2}_top_cand.txt | awk -v k="$k" -v gravel="$gravel" -v selected="$selected" -v freq="$freq" '{print $1"\t"$2"\t"$7"\t"$8"\t"NR"\t"gravel"\t"selected"\tJSFS\tresamp_"k"\t"freq"\t"$6}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topTen_mastertable.txt

if [[ ${gen1} == ${gen2} ]] || [[ ! -z "$5" ]]; then
    # GWAS results
    counter=0
    while read line; do
      counter=$(echo ${counter}+1 | bc -l)
      #echo $counter
      pos=$(echo ${line} | cut -d'_' -f2)
      snp_ld=$(zcat ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/nonWF_${samples}_${samples}_${gen1}_${gen2}_FST.out.gz | mawk -v chr="$chr" -v pos="$pos" '$1==chr && $2<=pos && max<$2{max=$2;s=$0}END{print s}' - | cut -f6 -)
      tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/gwas_${gen1}_${gen2}_top_cand.txt | grep -w ${line} - | mawk -v snp_ld="$snp_ld" -v k="$k" -v gravel="$gravel" -v selected="$selected" -v counter="$counter" -v freq="$freq" '{print $1"\t"$3"\t"$2"\t"$11"\t"counter"\t"gravel"\t"selected"\tGWAS\tresamp_"k"\t"freq"\t"snp_ld}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topTen_mastertable.txt
    done < <(tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/gwas_${gen1}_${gen2}_top_cand.txt | cut -f2 -)
fi

#iHS results
#counter=0
#while read line; do
 # counter=$(echo ${counter}+1 | bc -l)
 # pos=$(echo ${line} | cut -d'_' -f2)
 # snp_ld=$(zcat ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/nonWF_${samples}_${gen1}_${gen2}_FST.out.gz | mawk -v chr="$chr" -v pos="$pos" '$1==chr && $2<=pos && max<$2{max=$2;s=$0}END{print s}' - | cut -f6 -)
 # tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/ihs_after_${gen2}.top_cand.txt | grep -w ${line} - | mawk -v snp_ld="$snp_ld" -v k="$k" -v gravel="$gravel" -v selected="$selected" -v counter="$counter" -v freq="$freq" '{print $2"\t"$3"\t"$1"\t"$8"\t"counter"\t"gravel"\t"selected"\tiHS\tresamp_"k"\t"freq"\t"snp_ld}' >> ${plots}/${model}_n${sample}_topTen_mastertable.txt
#done < <(tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/ihs_after_${gen2}.top_cand.txt | cut -f1 -)

done
done

### This is based on Top Ten but can be adjusted to whatever "Top" user wants.
awk '{ if ($3==$7 && $5==1) {print $0 "\tOne"} else if ($3==$7 && $5<=5) {print $0 "\tFive"} else if ($3==$7 && $5<=10) {print $0 "\tTen"} else {print $0"\tNone"} }' ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt > ${plots}/${gen1}_${gen2}/${model}_n${sample}_tmp
mv ${plots}/${gen1}_${gen2}/${model}_n${sample}_tmp ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt
echo "Table complete."
