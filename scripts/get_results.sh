#!/bin/bash

echo -e "\nSOS: Get results for multiple simulations and resamplings.\n"
# usage-message
: ${1?"Usage: $0 <model> <generation1> <generation2> <sample_size> [demeGeneration] [midchrom] [nsims] [nresamp] [workdir]"
"NOTE: If different sample sizes per generation/group or multiple chromosomes, provide '_' between values "
"e.g. bash $0 add_grad_v1.0_f0.1_sim 40001 40003 10_15"
}

if [[ $# -lt 4 && $# -gt 0 ]] ; then
    echo -e 'You are missing arguments.\nUsage: bash power_calculate.sh <model> <generation1> <generation2> <sample_size> [demeGeneration] [midchrom] [nsims] [nresamp] [workdir]'
    exit 0
elif [[ $# -gt 4 && $# -lt 10 ]] ; then
    echo 'You can either provide all 10 arguments OR the first 4 with adjusted parameters within the script.'
    exit 0
elif [[ $# -eq 4 ]] ; then
    ## CHANGE THESE PARAMETERS FOR USER NEEDS
    deme=40000
    midchrom=33759515
    nsims=10
    nresamp=100
    workd="/home/lnd775/data/SimOutbreakSelection/rinderpest/outputs/"
    chrom="21_22"
else
    deme=$5
    midchrom=$6
    nsims=$7
    nresamp=$8        
    chrom=$9
    workd=${10}
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
workd=${workd}"/VAA"${VAA}/

if [[ ${sample} == *"_"* ]]; then
	samples=${sample}
	sample1=$(cut -d'_' -f1 <<< ${sample})
	sample2=$(cut -d'_' -f2 <<< ${sample})
	sample=$(echo ${sample1}"_n"${sample2})
else
	samples=${sample}
	sample=$(echo ${sample}"_n"${sample})
fi

if [ ! -d "$plots/" ]; then
	mkdir -p ${plots}/
fi

cd $plots/

if [ ! -d "${gen1}_${gen2}" ]; then
    mkdir ${gen1}_${gen2}
fi

for ((i=1;i<=$nsims;i++)); do
echo $i
sim=${i}
gravel=$sim
for ((k=1;k<=$nresamp;k++)); do
selected=$(sed -n 1p ${workd}${model}${i}/selected.txt)
freq=$(grep -w ${deme} ${workd}${model}${i}/${model}${i}.out | sed -n 2p | awk '{print $3}')

echo ${midchrom}
echo ${chrom}

if [[ ${midchrom} == *"_"*  || ${chrom} == *"_"* ]]; then
    echo "YUP"
    readarray -td '_' chrom_array < <(printf '%s' "$chrom");
    #unset 'chrom_array[-1]'; 
    declare -p chrom_array;
    readarray -td '_' midchrom_array < <(printf '%s' "$midchrom"); 
    #unset 'midchrom_array[-1]'; 
    declare -p midchrom_array;
else
    echo "NOPE"
	midchrom_array=(${midchrom})
    chrom_array=(${chrom})
fi

# Print arrays for debugging
echo "Chromosome Array: ${chrom_array[@]}"
echo "Midchrom Array: ${midchrom_array[@]}"
echo "Selected position: $selected"

# Loop over the midchrom and chrom arrays to determine which chromosome the selected position belongs to
for b in "${!midchrom_array[@]}"; do
    mid=${midchrom_array[$b]}
    next_mid=${midchrom_array[$((b+1))]:-9999999999}  # Handle the last chromosome by giving a large number

    # Debug prints for each step
    echo "Checking: Selected=$selected, Mid=$mid, Next Mid=$next_mid, Index=$b"

    if [[ $selected -ge $mid && $selected -lt $next_mid ]]; then
        # Subtract mid value and assign the corresponding chromosome
        selected=$(echo "$selected - $mid" | bc -l)
        selected="chr${chrom_array[$((b+1))]}_$selected"
        chr=${chrom_array[$((b+1))]}
        echo "Assigned to Chromosome: $selected"
        break
    elif [[ $selected -lt ${midchrom_array[0]} ]]; then
        # If selected position is less than the first midchrom, assign to first chromosome
        selected="chr${chrom_array[0]}_$selected"
        chr=${chrom_array[0]}
        echo "Assigned to First Chromosome: $selected"
        break
    fi
done

# FST results
tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/fst_${gen1}_${gen2}_top_cand.txt | awk -v k="$k" -v gravel="$gravel" -v selected="$selected" -v freq="$freq" '{print $1"\t"$2"\t"$7"\t"$3"\t"NR"\t"gravel"\t"selected"\tFST\tresamp_"k"\t"freq"\t"$6}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt
# JSFS results
tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/jsfs_${gen1}_${gen2}_top_cand.txt | awk -v k="$k" -v gravel="$gravel" -v selected="$selected" -v freq="$freq" '{print $1"\t"$2"\t"$7"\t"$8"\t"NR"\t"gravel"\t"selected"\tJSFS\tresamp_"k"\t"freq"\t"$6}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt

if [[ ${gen1} == ${gen2} ]]; then
    # GWAS results
    counter=0
    while read line; do
      counter=$(echo ${counter}+1 | bc -l)
      #echo $counter
      pos=$(echo ${line} | cut -d'_' -f2)
      snp_ld=$(zcat ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/nonWF_${samples}_${samples}_${gen1}_${gen2}_FST.out.gz | mawk -v chr="$chr" -v pos="$pos" '$1==chr && $2<=pos && max<$2{max=$2;s=$0}END{print s}' - | cut -f6 -)
      tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/gwas_${gen1}_${gen2}_top_cand.txt | grep -w ${line} - | mawk -v snp_ld="$snp_ld" -v k="$k" -v gravel="$gravel" -v selected="$selected" -v counter="$counter" -v freq="$freq" '{print $1"\t"$3"\t"$2"\t"$11"\t"counter"\t"gravel"\t"selected"\tGWAS\tresamp_"k"\t"freq"\t"snp_ld}' >> ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt
    done < <(tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/gwas_${gen1}_${gen2}_top_cand.txt | cut -f2 -)
fi

#iHS results
for name in ihs*; do
    if [[ -f "$name" ]]; then
        counter=0
        while read line; do
            counter=$(echo ${counter}+1 | bc -l)
            pos=$(echo ${line} | cut -d'_' -f2)
            snp_ld=$(zcat ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/nonWF_${samples}_${gen1}_${gen2}_FST.out.gz | mawk -v chr="$chr" -v pos="$pos" '$1==chr && $2<=pos && max<$2{max=$2;s=$0}END{print s}' - | cut -f6 -)
            tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/ihs_after_${gen2}.top_cand.txt | grep -w ${line} - | mawk -v snp_ld="$snp_ld" -v k="$k" -v gravel="$gravel" -v selected="$selected" -v counter="$counter" -v freq="$freq" '{print $2"\t"$3"\t"$1"\t"$8"\t"counter"\t"gravel"\t"selected"\tiHS\tresamp_"k"\t"freq"\t"snp_ld}' >> ${plots}/${model}_n${sample}_topCand_mastertable.txt
        done < <(tail -n +2 ${workd}${model}${i}/n${sample}/${k}/${gen1}_${gen2}/ihs_after_${gen2}.top_cand.txt | cut -f1 -)
    #else
    #    echo "No iHS files here."
    fi
done

done
done

### This is based on Top Ten but can be adjusted to whatever "Top" user wants.
awk '{ if ($3==$7 && $5==1) {print $0 "\tOne"} else if ($3==$7 && $5<=5) {print $0 "\tFive"} else if ($3==$7 && $5<=10) {print $0 "\tTen"} else {print $0"\tNone"} }' ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt > ${plots}/${gen1}_${gen2}/${model}_n${sample}_tmp
mv ${plots}/${gen1}_${gen2}/${model}_n${sample}_tmp ${plots}/${gen1}_${gen2}/${model}_n${sample}_topCand_mastertable.txt
echo "Table complete."
