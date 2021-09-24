#!/bin/bash

# Put in the start and end of each replicate
start=$1
end=$2

nonWF () {

#Set Work Directory
cd ~/results/selection/litterSize/
#cd ~/results/selection/

#Call SLiM simulation
slim -d I=${i} ~/epiSelect/simulations/originals/final_nonWFSelect.slim > rep_halffitness_${i}.out
#slim -d I=${i} ~/epiSelect/simulations/selection/final_nonWFSelectLitter.slim > rep_halffitness_${i}.out

# Takes nonWF SLiM outputs and transposes them into columns
for ((gen=58001;gen<=58009;gen++)); do
geng=${gen}g

# For nonWF outputs: putative parents with genotypes, mated parents, survivor parents
gawk -v gen="$gen" '$1==gen {print $0}' nwf_putparents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - > nwf_tmp${i}
gawk -v geng="$geng" '$1==geng {print $0}' nwf_putparents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' nwf_tmp${i} - | awk -v gen="$gen" '{print gen"\t"$1"\t"$2}' >> nwf_all_parents_${i}.out

gawk -v gen="$gen" '$1==gen {print $0}' nwf_mated_parents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - | awk -v gen="$gen" '{print gen"\t"$1}' >> nwf_mated_parents_${i}.out
gawk -v gen="$gen" '$1==gen {print $0}' nwf_survivor_parents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - | awk -v gen="$gen" '{print gen"\t"$1}' >> nwf_survivor_parents_${i}.out

gawk -v gen="$gen" '$1==gen {print $0}' nwf_mated_offspring_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - > nwf_tmp2_${i}
gawk -v geng="$geng" '$1==geng {print $0}' nwf_mated_offspring_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' nwf_tmp2_${i} - | awk -v gen="$gen" '{print gen"\t"$1"\t"$2}' >> nwf_mated_offspring_${i}.out

rm nwf_tmp${i}
rm nwf_tmp2_${i}
done
}

# Loop to create nonWF simulation replicates which produce selection-driven Nes
for ((i=${start};i<=${end};i++)); do
nonWF &
done
wait

WF () {
# Information regarding Ne from nonWF output are supplied to the WF models
cd ~/results/selection/litterSize/
#cd ~/results/selection/

z="$(tail -17 rep_halffitness_${i}.out | grep -w '58001' | cut -f2 | sed -n 1p)"

a="$(tail -17 rep_halffitness_${i}.out | grep -w '58002' | cut -f2 | sed -n 1p)"

b="$(tail -17 rep_halffitness_${i}.out | grep -w '58003' | cut -f2 | sed -n 1p)"

c="$(tail -17 rep_halffitness_${i}.out | grep -w '58004'| cut -f2 | sed -n 1p)"

d="$(tail -17 rep_halffitness_${i}.out | grep -w '58005'| cut -f2 | sed -n 1p)"

e="$(tail -17 rep_halffitness_${i}.out | grep -w '58006'| cut -f2 | sed -n 1p)"

f="$(tail -17 rep_halffitness_${i}.out | grep -w '58007'| cut -f2 | sed -n 1p)"

g="$(tail -17 rep_halffitness_${i}.out | grep -w '58008'| cut -f2 | sed -n 1p)"

h="$(tail -17 rep_halffitness_${i}.out | grep -w '58009'| cut -f2 | sed -n 1p)"

slim -d Z=${z} -d A=${a} -d B=${b} -d C=${c} -d D=${d} -d EE=${e} -d FF=${f} -d G=${g} -d H=${h} -d I=${i} ~/epiSelect/simulations/selection/paired/WFSelect_st_nwf_rep.slim > WF_112st_nwf_${i}.out

# Takes WF SLiM outputs and transposes them into columns
for ((gen=58001;gen<=58009;gen++)); do
geng=${gen}g

# For WF outputs: putative parents with genotypes and survivor parents
gawk -v gen="$gen" '$1==gen {print $0}' WF_putparents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - > wf_tmp${i}
gawk -v geng="$geng" '$1==geng {print $0}' WF_putparents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' wf_tmp${i} - | awk -v gen="$gen" '{print gen"\t"$1"\t"$2}' >> WF_all_parents_${i}.out

gawk -v gen="$gen" '$1==gen {print $0}' WF_result_parents_${i}.csv | cut -f2- | tr '\t' '\n' | paste -d'\t' - | awk -v gen="$gen" '{print gen"\t"$1}' >> WF_real_parents_${i}.out
rm wf_tmp${i}
done
}

# Loop which takes the Nes from the nonWF during the bottleneck and inserts them into a WF model equivalent and creates paired replicates
for ((i=${start}; i<=${end}; i++)); do
WF &
done
wait
