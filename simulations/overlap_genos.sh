#!/bin/bash

#Put in the start and end of each replicate
start=$1
end=$2

#cd ~/results/selection/
cd ~/results/selection/litterSize/
mv *.csv interim/
mv rep_halffitness_* slim_info/.
mv WF_112st_* slim_info/.


# nonWF data extraction
geno_cat () {

for ((gen=58001;gen<=58009;gen++)); do
# Finds the survivor parent IDs, extracts genotype, count number of offspring
cat nwf_survivor_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq - > nwf_tmp_${i}

while read line; do grep -h -w ${line} nwf_all_parents_${i}.out  >> nwf_survivor_patGeno${i}.log; done < <(cat nwf_tmp_${i} | cut -f2)

cat nwf_survivor_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq -c - >> nwf_survivor_patCount_${i}.log

rm  nwf_tmp_${i}

# Finds the mated parent IDs, extracts genotype, count number of offspring
cat nwf_mated_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq - > nwf_tmp_${i}

while read line; do grep -h -w ${line} nwf_all_parents_${i}.out  >> nwf_mated_patGeno${i}.log; done < <(cat nwf_tmp_${i} | cut -f2)

cat nwf_mated_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq -c - >> nwf_mated_patCount_${i}.log

rm  nwf_tmp_${i}

done
}


# WF data extraction
WF_geno_cat () {

for ((gen=58001;gen<=58009;gen++)); do

cat WF_real_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq - > wf_tmp_${i}

while read line; do grep -h -w ${line} WF_all_parents_${i}.out  >> WF_trupatGeno${i}.log; done < <(cat wf_tmp_${i} | cut -f2)

cat WF_real_parents_${i}.out | grep -w ${gen} | sort -k2 | uniq -c - >> WF_realCount_${i}.log

rm  wf_tmp_${i}

done
}


combi_counts_gt () {
#For nonWF outputs
#cd ~/results/selection/
cd ~/results/selection/litterSize/

awk '{print $1}' nwf_mated_patCount_${i}.log | paste nwf_mated_patGeno${i}.log - > processed/nwf_mated_${i}.log
awk '{print $1}' nwf_survivor_patCount_${i}.log | paste nwf_survivor_patGeno${i}.log - > processed/nwf_survivor_${i}.log

# Find IDs of individuals who did not mate - No reproduction
awk '{print $3}' nwf_mated_patCount_${i}.log | grep -v -f - nwf_all_parents_${i}.out > processed/nwf_nobabies_${i}.log

# Find parent IDs of those who died in the epidemic - Epidemic dead
awk '{print $3}' nwf_survivor_patCount_${i}.log | grep -v -f - processed/nwf_mated_${i}.log > processed/nwf_epidemic_${i}.log

#For WF outputs
awk '{print $1}' WF_realCount_${i}.log | paste WF_trupatGeno${i}.log - > processed/WF_mated_${i}.log

# No reproduction
awk '{print $3}' WF_realCount_${i}.log | grep -v -f - WF_all_parents_${i}.out > processed/WF_nobabies_${i}.log
}

concat_gts () {
#Concatenate all simulated iterations
cat processed/nwf_mated_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/nwf_mated.log
cat processed/nwf_survivor_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/nwf_survivor.log
cat processed/nwf_nobabies_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/nwf_nobabies.log
cat processed/nwf_epidemic_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/nwf_epidemic.log
cat processed/WF_mated_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/WF_mated.log
cat processed/WF_nobabies_${i}.log | awk -v i="$i" '{print $0"\t"i}' >> processed/WF_nobabies.log
mv processed/*_${i}.log processed/iter/.
}


for ((i=${start};i<=${end};i++)); do
geno_cat &
WF_geno_cat &
done
wait

for ((i=${start};i<=${end};i++)); do
combi_counts_gt
concat_gts
done
wait

mv *.out logouts/.
mv *patGeno*.log logouts/.
mv *Count*.log logouts/.
