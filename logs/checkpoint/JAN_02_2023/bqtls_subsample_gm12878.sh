dsqtl=results/variant_data/bqtls/pu1/pu1_meta_data.tsv
genome=reference/male.hg19.fa


model_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0
model=$model_dir/chrombpnet_model/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0



model_dir=results/chrombpnet/ATAC_PE/GM12878/subsampling/GM12878_250M/
model=$model_dir/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0


model_dir=results/chrombpnet/ATAC_PE/GM12878/subsampling/GM12878_100M/
model=$model_dir/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0


model_dir=results/chrombpnet/ATAC_PE/GM12878/subsampling/GM12878_50M/
model=$model_dir/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0



model_dir=results/chrombpnet/ATAC_PE/GM12878/subsampling/GM12878_25M/
model=$model_dir/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0

model_dir=results/chrombpnet/ATAC_PE/GM12878/subsampling/GM12878_5M/
model=$model_dir/chrombpnet_wo_bias.h5

output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
