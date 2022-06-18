chrombpnet_nb=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/chrombpnet_model/bias_model_scaled.h5
celline=K562
gpu=7
output_dir=results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/interpret/

#zcat $regions > $output_dir/peaks.bed
regions=imp_k562_new.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


#results/chrombpnet_feb_04/DNASE_SE/K562/DNASE_SE_01.24.2022_+0.9/chrombpnet_model/

