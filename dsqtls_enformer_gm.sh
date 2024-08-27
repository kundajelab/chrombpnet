dsqtl=/mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/dsqtl_meta_data.tsv
genome=/mnt/lab_data2/anusri/chrombpnet/reference/male.hg19.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
output_dirn=/mnt/lab_data2/anusri/variant-scorer/src/output/dsqtls_lcl/enformer_preds_small_window/
gpu=MIG-166d7783-762d-5f61-b31c-549eb4e0fba0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer_new_center.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0


