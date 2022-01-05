## MAKE FOOTPRINTS

cell_line=K562
data_type="ATAC"
neg_shift=4

date=$(date +'%m.%d.%Y')
setting=4_$neg_shift"_shifted_"$data_type"_"$date"_subsample_100M"
cur_file_name="withk562bias_k562_atac_subsample_100M.sh"
setting=4_4_shifted_ATAC_09.30.2021_subsample_100M
### SIGNAL INPUT

in_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/subsampled_ATAC_K562/bulk/K562.filtered.merged.1.82e-1.bam
overlap_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/K562/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optima$
idr_peak=/oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/K562/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narro$
is_filtered=True
samtools_flag=None

blacklist_region=$PWD/data/GRch38_unified_blacklist.bed
chrom_sizes=$PWD/data/hg38.chrom.sizes
ref_fasta=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

main_dir=$PWD/results/chrombpnet/$data_type/$cell_line
data_dir=$PWD/results/chrombpnet/$data_type/$cell_line/data_subsample_100M

### MODEL PARAMS

gpu=2
seed=1234
model_name=model
neg_dir=$main_dir/negatives_data
flank_size=1057
n_dil_layers=8
filters=500
neg_bed_test=$neg_dir"/bpnet.inputs.test.0.negatives.bed"

setting=4_4_shifted_ATAC_09.29.2021_bias_filters_500
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set1 --model_dir $output_dir/final_model_step3/unplug/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set2 --model_dir $output_dir/final_model_step3/unplug/


setting=4_4_shifted_ATAC_09.30.2021_subsample_100M
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set1 --model_dir $output_dir/with_k562_bias_final_model/unplug/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set2 --model_dir $output_dir/with_k562_bias_final_model/unplug/


setting=4_4_shifted_ATAC_10.01.2021_subsample_50M
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set1 --model_dir $output_dir/with_k562_bias_final_model/unplug/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set2 --model_dir $output_dir/with_k562_bias_final_model/unplug/

setting=4_4_shifted_ATAC_09.30.2021_subsample_25M
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set1 --model_dir $output_dir/with_k562_bias_final_model/unplug/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set2 --model_dir $output_dir/with_k562_bias_final_model/unplug/


setting=4_4_shifted_ATAC_10.01.2021_subsample_5M
output_dir=$PWD/results/chrombpnet/$data_type/$cell_line/$setting

CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set1 --model_dir $output_dir/with_k562_bias_final_model/unplug/
CUDA_VISIBLE_DEVICES=$gpu python  $PWD/src/evaluation/marginal_footrprints/main_footprints_check.py --ref_fasta $ref_fasta --gc_neg $neg_bed_test --motifs k562_motifs_set2 --model_dir $output_dir/with_k562_bias_final_model/unplug/


