dsqtl=/mnt/lab_data2/anusri/chrombpnet/results/variant_data/dsqtls/dsqtl_meta_data.tsv
#dsqtl=/mnt/lab_data2/anusri/chrombpnet/results/variant_data/dsqtls/test.tsv

genome=/mnt/lab_data3/anusri/histone_expts/all_qtl_analysis/dnase_qtl/data/hg18.fa
#model=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
model=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5

#output_dirn=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/dsqtl_preds/
output_dirn=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/dsqtl_preds_enformer_test/
mkdir $output_dirn
gpu=0

#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction/snp_scoring_enformer.py -i $dsqtl -g $genome  -o $output_dirn -bs 1 --debug_mode_on 0


