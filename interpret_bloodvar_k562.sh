dsqtl=/mnt/lab_data2/anusri/signed_variant_scorer/variant-scorer/output/blood_traits/interpret_bloodtraits.tsv

#blood_trait_variants.tsv
genome=reference/hg38.genome.fa

model_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/
model=$model_dir/chrombpnet_model/chrombpnet_wo_bias.h5


output_dirn=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/blood_traits_preds/
mkdir $output_dirn

gpu=2

#CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
#mkdir $output_dirn/bigwigs/
#python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes reference/chrom.sizes 


genome=reference/hg38.genome.fa

model_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/
model=$model_dir/chrombpnet_model/chrombpnet_wo_bias.h5


output_dirn=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/blood_traits_preds_dnase/
mkdir $output_dirn

gpu=2

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
mkdir $output_dirn/bigwigs/
python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes reference/chrom.sizes 



