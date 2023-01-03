dsqtl=results/variant_data/bqtls/pu1/pu1_meta_data.tsv
genome=reference/male.hg19.fa


model_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0
model=$model_dir/chrombpnet_model/chrombpnet_wo_bias.h5


output_dirn=$model_dir/bqtls_pu1_preds/
mkdir $output_dirn

gpu=0

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0


