chrombpnet_nb=/oak/stanford/groups/akundaje/anusri/Mineto-data/MemoryCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=/oak/stanford/groups/akundaje/anusri/Mineto-data/MemoryCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/chrombpnet.h5
bias=/oak/stanford/groups/akundaje/anusri/Mineto-data/MemoryCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/bias_model_scaled.h5
celline=MemoryCD4
gpu=2

#main_dir=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/
#output_dir=$main_dir/interpret/

output_dir=/oak/stanford/groups/akundaje/anusri/Mineto-data/MemoryCD4/ATAC_PE_09.21.2022_withgm12878bias/chrombpnet_model/interpret_full/
mkdir $output_dir

regions=/oak/stanford/groups/akundaje/anusri/Mineto-data/MemoryCD4/data/peaks_no_blacklist.bed 
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


