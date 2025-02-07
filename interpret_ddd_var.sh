dsqtl=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/glutameric_neuron/snps.bed


genome=reference/hg38.genome.fa

model=/oak/stanford/groups/akundaje/salil512/chd/model_training/chrombpnet/fetal_brain/chrombpnet_models/c1/fold_0/chrombpnet_wo_bias.h5



output_dirn=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/glutameric_neuron/fold0/preds/
mkdir $output_dirn


gpu=2

CUDA_VISIBLE_DEVICES=$gpu python src/evaluation/variant_effect_prediction_interpret/snp_scoring.py -i $dsqtl -g $genome -m $model -o $output_dirn -bs 64 --debug_mode_on 0
mkdir $output_dirn/bigwigs/
python src/evaluation/variant_effect_prediction_interpret/convert_shap_to_bigiwg.py  -p $output_dirn -o $output_dirn/bigwigs/ --chromsizes reference/chrom.sizes 



