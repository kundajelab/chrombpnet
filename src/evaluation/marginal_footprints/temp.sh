cell_type="GM12878"

model="/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0"

output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=0 python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        -r $regions \
        -fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
        -cm $model_dir"/chrombpnet_wo_bias.h5" \
        -um $model_dir"/chrombpnet.h5" \
        -bs 128 \
        -o $output$cell_type \
        -pwm_f union_motifs.tsv


