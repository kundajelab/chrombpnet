chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
bias=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=1
output_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/bias/

#merge gm12878 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/optimal_overlap_peaks/GM12878.overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | uniq > results/chrombpnet/gm12878.merged.atac.dnase.peaks.bed

#zcat $regions > $output_dir/peaks.bed
#regions=gm_locus_peaks.bed
regions=combine_gm.bed
bash make_bigwig_new_w_bias.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


