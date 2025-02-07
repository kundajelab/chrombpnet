
mkdir output/ATAC
mkdir output/DNASE_NEW

cell_type=IMR90
dtyle=DNASE_SE

mkdir output/ATAC/$cell_type
mkdir output/DNASE_NEW/$cell_type

CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds_new.py \
        -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        -c $cell_type \
        -o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
        -r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
        -pwm_f union_motifs_fig.tsv


cell_type=GM12878
dtyle=ATAC_PE

mkdir output/ATAC/$cell_type
mkdir output/DNASE_NEW/$cell_type

CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds_new.py \
	-g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-c $cell_type \
	-o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
	-pwm_f union_motifs_fig.tsv

cell_type=H1ESC
dtyle=ATAC_PE

mkdir output/ATAC/$cell_type
mkdir output/DNASE_NEW/$cell_type

CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds_new.py \
	-g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
	-c $cell_type \
	-o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
	-pwm_f union_motifs_fig.tsv




cell_type=HEPG2
dtyle=ATAC_PE

mkdir output/ATAC/$cell_type
mkdir output/DNASE_NEW/$cell_type


#CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds_new.py \
#	-g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
#	-c $cell_type \
#	-o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
#	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
#	-pwm_f union_motifs_fig.tsv

cell_type=K562
dtyle=ATAC_PE

mkdir output/ATAC/$cell_type
mkdir output/DNASE_NEW/$cell_type

#CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds_new.py \
#	-g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
#	-c $cell_type \
#	-o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
#	-r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
#	-pwm_f union_motifs_fig.tsv

