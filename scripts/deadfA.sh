workd="/home/lnd775/data/plague/results/sweden/VAA08/"
#model="add_grad_s0.1_f0.2_sim"

for model in rec_grad_s0.1_f0.1_sim rec_grad_s0.1_f0.2_sim rec_grad_s0.1_f0.3_sim rec_grad_s0.1_f0.4_sim; do
for i in {1..10}; do

cd ${workd}${model}${i}/
selected=$(sed -n 1p ${workd}${model}${i}/selected.txt)

python3 /home/lnd775/data/plague/scripts/just_dead_af.py --tree ${workd}${model}${i}/rec_grad_${i}.trees --muTree ${workd}${model}${i}/rec_grad_${i}_mu.trees --bottle 58799 --dead 58799 --select ${selected} --generations 58799 58799 --putparents ${workd}${model}${i}/nwf_putparents_${i}.csv --sampling 10 10 --nosel --deme 58796 --seed 101 > ${workd}${model}${i}/dead_af_${i}.txt

dead=$(grep -A1 'DEAD' dead_af_${i}.txt | tail -n1 | awk '{print $1}')
python3 /home/lnd775/data/plague/scripts/just_dead_af.py --tree ${workd}${model}${i}/rec_grad_${i}.trees --muTree ${workd}${model}${i}/rec_grad_${i}_mu.trees --bottle 58799 --dead 58799 --select ${selected} --generations 58799 58799 --putparents ${workd}${model}${i}/nwf_putparents_${i}.csv --sampling ${dead} 10 --nosel --deme 58796 --seed ${i}

if [[ $selected -ge 33759515 ]]; then
  selected=$(echo $selected-33759515 | bc -l)
  selected="chr22_"$selected
else
  selected="chr21_"$selected
fi
echo $selected;

zgrep -m 1 -w ${selected} nonWF_${dead}_10_58799_58799_FST.out.gz | cut -f4 | awk -v dead="$dead" '{print $1/(dead*2)}' - >> plague_dead_readAF.txt
done
done

