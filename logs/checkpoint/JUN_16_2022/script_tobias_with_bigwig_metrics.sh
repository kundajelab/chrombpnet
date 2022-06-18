python /mnt/lab_data2/anusri/chrombpnet/src/training/predict.py \
        --genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        --bigwig=results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw \
        --peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed \
        --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
        --inputlen=2114 \
        --outputlen=1000 \
        --output_prefix=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias_wo_bias \
        --batch_size=256 \
        --model_h5=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias_wo_bias.h5

python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_with_tobias.py \
        --genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        --bigwig=results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw \
        --bias_bigwig=results/tobias/ATAC_PE/GM12878/data/GM12878_expected_not_multiplied_with_counts_2/sorted_merged_expected.bw \
        --peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed \
        --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
        --inputlen=2114 \
        --outputlen=1000 \
        --output_prefix=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias \
        --batch_size=256 \
        --model_h5=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias.h5


python /mnt/lab_data2/anusri/chrombpnet/src/training/predict_with_bigwig.py \
        --genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa \
        --bigwig=results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw \
        --bias_bigwig=results/tobias/ATAC_PE/GM12878/data/GM12878_expected_not_multiplied_with_counts_2/sorted_merged_expected.bw \
        --peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed \
        --nonpeaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed \
        --chr_fold_path=/mnt/lab_data2/anusri/chrombpnet/splits/fold_0.json \
        --inputlen=2114 \
        --outputlen=1000 \
        --output_prefix=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_bias_model/tobias \
        --batch_size=256 \
        --model_h5=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_bias_model/tobias_bias.h5
