
#split -l 78366 /mnt/lab_data2/anusri/variant-scorer/src/output/caqtls_lcl_latest/enformer_preds_small_window/metad_data.tsv /mnt/lab_data2/anusri/variant-scorer/src/output/caqtls_lcl_latest/enformer_preds_small_window/split

dsqtl=/mnt/lab_data2/anusri/variant-scorer/src/output/caqtls_lcl_latest/enformer_preds_small_window/splitab
genome=/mnt/lab_data2/anusri/chrombpnet/reference/male.hg19.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
output_dirn=/mnt/lab_data2/anusri/variant-scorer/src/output/caqtls_lcl_latest/enformer_preds_small_window/splitab/
mkdir $output_dirn

gpu=2

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer_new_center.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0



