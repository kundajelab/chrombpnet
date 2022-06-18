chrombpnet_nb=results/chrombpnet_feb_04/ATAC_PE/HEPG2/ATAC_PE_withk562bias_12.30.2021/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet_feb_04/ATAC_PE/HEPG2/ATAC_PE_withk562bias_12.30.2021/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet_feb_04/ATAC_PE/HEPG2/ATAC_PE_withk562bias_12.30.2021/chrombpnet_model/bias_model_scaled.h5
celline=HEPG2
gpu=0
regions=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/HEPG2/peaks.bed.gz
output_dir=results/chrombpnet_feb_04/ATAC_PE/HEPG2/ATAC_PE_withk562bias_12.30.2021/interpret/

#zcat $regions > $output_dir/peaks.bed
regions=$output_dir/peaks.bed
bash make_bigwig.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


