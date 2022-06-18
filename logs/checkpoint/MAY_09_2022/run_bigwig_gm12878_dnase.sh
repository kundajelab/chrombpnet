chrombpnet_nb=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5 
chrombpnet=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=0
output_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/interpret/

#merge gm12878 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/optimal_overlap_peaks/GM12878.overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | uniq > results/chrombpnet/gm12878.merged.atac.dnase.peaks.bed

#zcat $regions > $output_dir/peaks.bed
regions=results/chrombpnet/gm12878.merged.atac.dnase.peaks.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


