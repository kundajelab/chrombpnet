chrombpnet_nb=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/chrombpnet_model/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/chrombpnet_model/bias_model_scaled.h5
celline=THP1
gpu=0
output_dir=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/chrombpnet_model/locus_interpret/


#zcat $regions > $output_dir/peaks.bed
regions=/oak/stanford/groups/akundaje/anusri/THP1-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_03.24.2022_withgm12878bias/chrombpnet_model/locus_interpret/peaks_in_main_locus.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


