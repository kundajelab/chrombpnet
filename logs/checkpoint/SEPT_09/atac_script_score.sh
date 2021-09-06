#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=H1
data_type="ATAC"

date=$(date +'%m.%d.%Y')
setting=$data_type"_07.22.2021"
cur_file_name="h1_atac_script.sh"

### SIGNAL INPUT


blacklist_region=$PWD/data/all_three_blacklists.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/$cell_line
data_dir=$PWD/$cell_line/data
output_dir=$PWD/$cell_line/$setting


### MODEL PARAMS

gpu=0
filters=500 
n_dil_layers=8
seed=1234 
model_name=model 
neg_dir=$main_dir/neg_data
flank_size=1057


fold=0
./main_scripts/invivo_bias_model_step1/score.sh $output_dir/invivo_bias_model_step1 $model_name $fold $cell_line $seed
./main_scripts/bias_fit_on_signal_step2/score.sh $output_dir/bias_fit_on_signal_step2 $model_name $fold $cell_line $seed


