
get_trufA() {
workd="/home/lnd775/data/plague/results/sweden/VAA${VAA}/"

# This is for the True allele frequency
freq=$(grep -w ${deme} ${workd}${model}${i}/${model}${i}.out | tail -n 1 | awk '{print $3}')
if [[ $gen1 == $gen2 ]]; then
    grep -w ^${gen1} ${workd}${model}${i}/${model}${i}.out | tail -n 1 | cut -f1-3 | awk -v rec="$rec" -v i="$i" '{print$0"\tsim_"i"\t"rec"\tFull"}' - >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
    cat ${workd}${model}${i}/plague_dead_readAF.txt | head -n 1 | cut -f1 | awk -v gen="$gen1" -v rec="$rec" -v i="$i" '{print gen"d\t0\t"$1"\tsim_"i"\t"rec"\tFull"}' - >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt

else
    grep -w ^${gen1} ${workd}${model}${i}/${model}${i}.out | tail -n 1 | cut -f1-3 | awk -v rec="$rec" -v i="$i" '{print$0"\tsim_"i"\t"rec"\tFull"}' - >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
    grep -w ^${gen2} ${workd}${model}${i}/${model}${i}.out | tail -n 1 | cut -f1-3 | awk -v rec="$rec" -v i="$i" '{print$0"\tsim_"i"\t"rec"\tFull"}' - >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
fi

}


# This is for estimated AF
get_freqs() {
workd="/home/lnd775/data/plague/results/sweden/VAA${VAA}/"

# Get the selected variant and assess which chromosome it's on
selected=$(sed -n 1p ${workd}${model}${i}/selected.txt)
if [[ ${selected} -ge ${midchrom} ]]; then
    chrom=22
    selected=$(echo ${selected}-${midchrom} | bc -l)
else
    chrom=21
fi

for ((k=1;k<=100;k++)); do
echo This is sim $i and resampling $k
if [[ $gen1 == $gen2 ]]; then
#Get counts for both Dead and Alive
	counts=$(zgrep -w chr${chrom}_${selected} ${workd}${model}${i}/n${sample}_n${sample}/${k}/${gen1}_${gen2}/nonWF_${sample}_${sample}_${gen1}_${gen2}_FST.out.gz | cut -f4-5)
	#echo ${counts}
	ac=$(echo -e ${counts} | awk '{print $1}')
    af=$(echo -e ${ac} | awk -v sample="$sample" '{print $0/(2*sample)}')
	#echo "Dead: "${ac}" and " ${af} 
    echo -e ${gen}"d\t"${sample}"\t"$af"\tsim_"$i"\t"${rec}"\tresamp_"${k} >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
    ac=$(echo -e ${counts} | awk '{print $2}')
    af=$(echo -e ${ac} | awk -v sample="$sample" '{print $0/(2*sample)}')
    #echo "Alive: "${ac}" and " ${af} 
    echo -e ${gen}"a\t"${sample}"\t"$af"\tsim_"$i"\t"${rec}"\tresamp_"${k} >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
else
	counts=$(zgrep -w chr${chrom}_${selected} ${workd}${model}${i}/n${sample}_n${sample}/${k}/${gen1}_${gen2}/nonWF_${sample}_${sample}_${gen1}_${gen2}_FST.out.gz | cut -f4-5)
	ac=$(echo -e ${counts} | awk '{print $1}')
	af=$(echo -e ${ac} | awk -v sample="$sample" '{print $0/(2*sample)}')
	echo -e ${gen1}"\t"${sample}"\t"$af"\tsim_"$i"\t"${rec}"\tresamp_"${k} >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt

   	ac=$(echo -e ${counts} | awk '{print $2}')
    af=$(echo -e ${ac} | awk -v sample="$sample" '{print $0/(2*sample)}')
    echo -e ${gen2}"\t"${sample}"\t"$af"\tsim_"$i"\t"${rec}"\tresamp_"${k} >> ${plots}/${gen1}_${gen2}/${model}_freq_n${sample}.txt
fi
done; }


for VAA in 08; do
for model in add_grad_s0.1_f0.4_sim; do
#for model in add_grad_s0.1_f0.1_sim add_grad_s0.1_f0.2_sim add_grad_s0.1_f0.3_sim add_grad_s0.1_f0.4_sim; do
#model="add_grad_s0.1_f0.1_sim"
ih=$(echo ${model} | cut -d'_' -f1)
midchrom=33759515
gen1=$1
gen2=$2
sample=$3
deme=58796
gen=$gen1
plots="/home/lnd775/data/plague/results/plots/sweden/${ih}/VAA${VAA}/"

cd $plots/
if [ ! -d "${gen1}_${gen2}" ]; then
    mkdir ${gen1}_${gen2}
fi

if [[ "$model" == *"grad"* ]]; then
    rec="Gradual"
#    end=17
else
    rec="Instant"
#    end=15
fi

for ((i=1;i<=10;i++)); do
get_freqs &
get_trufA &
echo "Freqs for Sim ${i} done."
done  

done &
done



