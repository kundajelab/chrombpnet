
gpu=3


cell_type="K562"
model="/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
        -tt $cell_type"_DNASE"



cell_type="H1ESC"
model="/nautilus_runs_apr12/H1ESC_04.09.2022_bias_128_4_1234_0.8_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
        -tt $cell_type"_DNASE"

cell_type="IMR90"
model="/nautilus_runs_apr12/IMR90_04.09.2022_bias_128_4_1234_0.8_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
        -tt $cell_type"_DNASE"





cell_type="HEPG2"


model="/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
        -tt $cell_type"_DNASE"


cell_type="GM12878"

model="/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0"
output="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/DNASE/"$cell_type"/footprints/"
regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/$cell_type/negatives_data/negatives_with_summit.bed
model_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/"$cell_type$model"/chrombpnet_model"
mkdir $output

CUDA_VISIBLE_DEVICES=$gpu python marginal_footprinting_with_uncorrected.py -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-r $regions \
	-fl /mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
	-cm $model_dir"/chrombpnet_wo_bias.h5" \
	-um $model_dir"/chrombpnet.h5" \
	-bs 128 \
	-o $output$cell_type \
	-pwm_f union_motifs.tsv \
        -tt $cell_type"_DNASE"
