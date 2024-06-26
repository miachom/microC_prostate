#!/bin/bash
source /etc/profile
#$ -S /bin/bash
#$ -pe pvm 8

source /home/ma1111/miniconda3/etc/profile.d/conda.sh
#source /home/ma1111/miniconda3/bin/activate
conda activate
conda init bash
conda activate /home/ma1111/miniconda3/envs/chromo_env


## in cchromisght example -m8000 -M50000 detect 8-80kb loop size with 
out_path=/home/ma1111/microC_project/prostate/compt_analysis/chromosight_run/loop
seq_path=/home/ma1111/microC_project/prostate/hic_nf_cools/cool_files
sample_list=('CT20_MicroC_1_S1' , 'CT20_MicroC_2_S2' , 'CT20_MicroC_3_S3' , 'CT20_MicroC_4_S4' , 'CT20_MicroC_5_S5' , 'CT20_MicroC_6_S6' , 'CT20_MicroC_7_S7' , 'CT20_MicroC_8_S8')

for i in "${sample_list[@]}"
do
for j in $seq_path/*.cool
do

chromosight detect --threads 8 ${j} -m100000 -M500000000 -p0.35 $out_path/${i}_50mb_loops
chromosight detect --threads 8 --pattern=borders --pearson=0.4 ${j} $out_path/${i}.50mb_borders
chromosight detect --threads 8 --pattern=hairpins --pearson=0.4 ${j} $out_path/${i}.50mb_hairpins

done
done

