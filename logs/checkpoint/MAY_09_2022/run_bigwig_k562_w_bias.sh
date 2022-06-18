chrombpnet_nb=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.5/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.5/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.5/chrombpnet_model/bias_model_scaled.h5
output_dir=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.5/interpret/
#results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.5/
celline=K562
gpu=2


regions=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/interpret/K562.interpreted_regions.bed
bash make_bigwig_new_w_bias.sh  $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


