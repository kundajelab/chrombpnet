#!/bin/bash
##TODO flank_size and ref_fasta files in scripts 

cell_line=K562
data_type="ATAC"
neg_shift=4
filters=500

date=$(date +'%m.%d.%Y')
setting=4_$neg_shift"_shifted_"$data_type"_"$date"_bias_filters_"$filters"_subsample_25M"
cur_file_name="k562_atac_subsample_25M.sh"
setting=4_4_shifted_ATAC_10.04.2021_bias_filters_500_subsample_25M
### SIGNAL INPUT
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/K562/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/K562/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
is_filtered=True
samtools_flag=None

blacklist_region=$PWD/data/GRch38_unified_blacklist.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$PWD/results/chrombpnet/$data_type/$cell_line/data_subsample_25M
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

### MODEL PARAMS

gpu=1
n_dil_layers=8
seed=1234 
model_name=model 
neg_dir=$main_dir/negatives_data
flank_size=1057


### STEP 1 -  TRAIN BIAS MODEL

min_logcount=0.0
max_logcount=2.5
model_name=model
fold=0
bash $PWD/src/models/chrombpnet_scripts/invivo_bias_model_step1/score.sh $output_dir/invivo_bias_model_with_500M_preds $model_name $fold $cell_line $seed $min_logcount $max_logcount

