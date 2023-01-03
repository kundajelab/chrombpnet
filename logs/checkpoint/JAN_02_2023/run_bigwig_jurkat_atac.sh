chrombpnet_nb=/oak/stanford/groups/akundaje/anusri/Jurkat-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_08.02.2022_withgm12878bias/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=/oak/stanford/groups/akundaje/anusri/Jurkat-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_08.02.2022_withgm12878bias/chrombpnet_model/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/anusri/Jurkat-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_08.02.2022_withgm12878bias/chrombpnet_model/bias_model_scaled.h5
celline=Jurkat
gpu=0


output_dir=/oak/stanford/groups/akundaje/anusri/Jurkat-Engreitz/results/chrombpnet/ATAC_PE/ATAC_PE_08.02.2022_withgm12878bias/interpret/
mkdir $output_dir


regions=jurkat_interest_peaks.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


