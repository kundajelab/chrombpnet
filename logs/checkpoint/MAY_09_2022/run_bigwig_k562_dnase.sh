chrombpnet_nb=results/chrombpnet/DNASE_PE/K562/nautilus_runs_apr12/K562_04.09.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet/DNASE_PE/K562/nautilus_runs_apr12/K562_04.09.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/DNASE_PE/K562/nautilus_runs_apr12/K562_04.09.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=K562
gpu=1
output_dir=results/chrombpnet/DNASE_PE/K562/nautilus_runs_apr12/K562_04.09.2022_bias_128_4_1234_0.8_fold_0/interpret/

#merge k562 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/optimal_overlap_peaks/K562.overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | uniq > results/chrombpnet/k562.merged.atac.dnase.peaks.bed
regions=results/chrombpnet/k562.merged.atac.dnase.peaks.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


