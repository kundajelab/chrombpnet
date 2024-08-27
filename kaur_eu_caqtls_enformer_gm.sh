
#split -l 33000 /mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/figure_6/final_figures/kumsaka.2018.kaur.scores.peak.filtered.bed  /mnt/lab_data2/anusri/variant-scorer/src/output/kaur_caqtls_lcl_latest/enformer_preds_small_window/splitss/split

dsqtl=/mnt/lab_data2/anusri/variant-scorer/src/output/kaur_caqtls_lcl_latest/enformer_preds_small_window/splitss/splitaa
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
output_dirn=/mnt/lab_data2/anusri/variant-scorer/src/output/kaur_caqtls_lcl_latest/enformer_preds_small_window/splitss/splitaa
mkdir $output_dirn

gpu=MIG-40f43250-998e-586a-ac37-d6520e92590f

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer_new_center.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0




