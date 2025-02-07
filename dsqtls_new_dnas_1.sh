#dsqtl=dsqtl_meta_data_small.tsv
dsqtl=dsqtl_meta_data_small_new.tsv
genome=reference/male.hg19.fa

model_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EMT/chrombpnet_model_feb15_fold_0/
model=$model_dir/chrombpnet_wo_bias.h5
output_dirn=$model_dir/dsqtls_interpret_small_1/
mkdir $output_dirn
gpu=3
chrom=reference/chrom.sizes 
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
mkdir $output_dirn/bigwigs/
#python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes $chrom



model_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/
model=$model_dir/chrombpnet_wo_bias.h5
output_dirn=$model_dir/dsqtls_interpret_small_new_1/
mkdir $output_dirn
gpu=3
chrom=reference/chrom.sizes 
#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
mkdir $output_dirn/bigwigs/
#python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes $chrom



output_dirn=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/over_fitting_test/enformer_dsqtls_interpret_small_new_1/
mkdir $output_dirn
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0


