
#tail -n +2 /mnt/lab_data2/anusri/enformer/eu_caqtls/source2.tsv > /mnt/lab_data2/anusri/variant-scorer/src/output/afr_caqtls_window/meta_data.tsv
#split -l 38833 /mnt/lab_data2/anusri/variant-scorer/src/output/afr_caqtls_window/meta_data.tsv  /mnt/lab_data2/anusri/variant-scorer/src/output/afr_caqtls_window/splits/split

dsqtl=/mnt/lab_data2/anusri/variant-scorer/src/output/afr_caqtls_window/splits/splitaa
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
#chrom_sizes=/mnt/data/annotations/by_release/hg19/hg19.chrom.sizes
output_dirn=/mnt/lab_data2/anusri/variant-scorer/src/output/afr_caqtls_window/splitaa/
mkdir $output_dirn
gpu=0


CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer_new_center.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0


