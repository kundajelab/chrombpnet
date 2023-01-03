chrombpnet_nb=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR020LUD/chrombppnet_model_encsr283tme_bias//chrombpnet_wo_bias.h5 
chrombpnet=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR020LUD/chrombppnet_model_encsr283tme_bias//chrombpnet.h5
bias=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR020LUD/chrombppnet_model_encsr283tme_bias//bias_model_scaled.h5
#chrombpnet_nb=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5 
#chrombpnet=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet.h5
#bias=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/bias_model_scaled.h5

celline=ENCSR020LUD
gpu=0
#main_dir=results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/
#output_dir=$main_dir/interpret/

output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/script/immune_scale/missing/
mkdir $output_dir

#merge k562 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/caper/K562.ENFF205FNC/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | uniq > results/chrombpnet/k562.merged.atac.dnase.peaks.bed

#regions=results/chrombpnet/k562.merged.atac.dnase.peaks.bed
regions=temp1.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


