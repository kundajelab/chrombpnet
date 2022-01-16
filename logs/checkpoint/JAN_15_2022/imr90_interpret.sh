
CUDA_VISIBLE_DEVICES=5 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/interpret/interpret.py  \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/IMR90/data/30K.subsample.overlap.bed    \
       --output_prefix=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/IMR90/ATAC_PE_12.30.2021/bias_model/interpret/IMR90    \
       --model_h5=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/IMR90/ATAC_PE_12.30.2021/bias_model/bias.h5 \
       --profile_or_counts=profile
