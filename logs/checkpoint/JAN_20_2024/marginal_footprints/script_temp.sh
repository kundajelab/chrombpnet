cell_type=IMR90
dtyle=DNASE_SE
CUDA_VISIBLE_DEVICES=1 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_with_uncorrected_all_folds.py \
        -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        -c $cell_type \
        -o /mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_20_2024/marginal_footprints/output \
        -r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$dtyle/$cell_type/negatives_data/negatives_with_summit.bed \
        -pwm_f union_motifs_fig.tsv
