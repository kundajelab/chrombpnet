chrombpnet_nb=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_2356/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_2356/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_2356/chrombpnet_model/bias_model_scaled.h5
celline=K562
gpu=0
regions=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz
output_dir=results/chrombpnet/ATAC_PE/K562_stability/K562_02.08.2022_bias_128_4_2356/interpret/

#zcat $regions > $output_dir/peaks.bed
regions=$output_dir/peaks.bed
bash make_bigwig.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


