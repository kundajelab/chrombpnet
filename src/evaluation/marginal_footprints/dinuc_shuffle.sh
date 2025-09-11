#CUDA_VISIBLE_DEVICES=2 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_all_dinuc_shuffle.py \
#        -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
#        -r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model//filtered.nonpeaks.bed \
#        --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
#        -m /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model//chrombpnet_wo_bias.h5 \
#        -bs 256 \
#        -o /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/all_motifs_footprints_shuffled/motif \
#        -pwm_f /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/gm_benchmarking_motifs.tsv \
#        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1


CUDA_VISIBLE_DEVICES=2 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/marginal_footprinting_all_dinuc_shuffle.py \
        -g /mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        -r /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model//filtered.nonpeaks.bed \
        --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
        -m /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model//chrombpnet_wo_bias.h5 \
        -bs 256 \
        -o /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/all_motifs_footprints_fwd/motif \
        -pwm_f /mnt/lab_data2/anusri/chrombpnet/src/evaluation/marginal_footprints/gm_benchmarking_motifs_new.tsv \
        -mo tn5_1,tn5_2,tn5_3,tn5_4,tn5_GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1

