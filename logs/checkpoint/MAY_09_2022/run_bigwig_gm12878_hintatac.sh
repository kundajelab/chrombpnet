chrombpnet_nb=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/hint_atac.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=2

output_dir=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/interpret/
#ATAC_PE_04.06.2022_tobias_with_bias_bigwig
#ATAC_PE_04.08.2022_tobias_corrected_not_softmax_custom_shift
#ATAC_PE_04.07.2022_tobias_corrected_not_softmax
#merge gm12878 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/GM12878/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/optimal_overlap_peaks/GM12878.overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | uniq > results/chrombpnet/gm12878.merged.atac.dnase.peaks.bed

#zcat $regions > $output_dir/peaks.bed
#regions=results/chrombpnet/gm12878.merged.atac.dnase.peaks.bed
regions=imp_locus_gm.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir


