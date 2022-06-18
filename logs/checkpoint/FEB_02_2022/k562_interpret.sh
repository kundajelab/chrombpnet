#regions=results/chrombpnet/ATAC_PE/K562/data/30K.subsample.overlap.bed
#regions=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/testing/beta_globin.coords
regions=/srv/scratch/anusri/colloboration_data_overlap/interesting.locus

mkdir /mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022/chrombpnet_model/interpret_no_weight
output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022/chrombpnet_model/interpret_no_weight
CUDA_VISIBLE_DEVICES=7 python /mnt/lab_data2/anusri/chrombpnet/src/evaluation/interpret/interpret_new.py  \
       --genome=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
       --regions=$regions  \
       --output_prefix=$output_dir/K562    \
       --model_h5=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/ATAC_PE_01.24.2022/chrombpnet_model/chrombpnet_wo_bias.h5 \
