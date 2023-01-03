chrombpnet_nb=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/chrombpnet_model/chrombpnet_wo_bias.h5
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=1


output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_12.03.2022_1234_8_2114_0_hepg2_transfer_bias/bias_interpret_gm_regions/
mkdir $output_dir

regions=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/30K.subsample.overlap.bed

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


