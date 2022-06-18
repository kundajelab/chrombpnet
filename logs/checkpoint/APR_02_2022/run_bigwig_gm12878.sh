chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=0
regions=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/peaks.bed.gz
output_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/interpret/

#zcat $regions > $output_dir/peaks.bed
regions=$output_dir/peaks.bed
bash make_bigwig.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


