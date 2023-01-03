dsqtl=/mnt/lab_data2/anusri/chrombpnet/results/variant_data/dsqtls/dsqtl_meta_data_small.tsv
genome=reference/hg18.fa

model_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0
model=$model_dir/chrombpnet_model/chrombpnet_wo_bias.h5


output_dirn=$model_dir/dsqtls_interpret_small/
mkdir $output_dirn

gpu=0

#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0


mkdir $output_dirn/bigwigs/
python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes /mnt/lab_data2/anusri/chrombpnet/results/variant_data/dsqtls/process/hg18.chrom.sizes

