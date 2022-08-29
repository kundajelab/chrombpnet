chrombpnet_nb=results/chrombpnet/ATAC_PE/K562/uncorrected_model_08.22.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
chrombpnet=results/chrombpnet/ATAC_PE/K562/uncorrected_model_08.22.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
bias=results/chrombpnet/ATAC_PE/K562/uncorrected_model_08.22.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
celline=K562
gpu=0

output_dir=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/uncorrected_model_08.22.2022_filters_512_dil_8/
mkdir $output_dir

#merge k562 peaks form atac and dnase
atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz
dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/caper/K562.ENFF205FNC/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
zcat $atac_peaks $dnase_peaks | bedtools sort -i stdin | uniq > results/chrombpnet/k562.merged.atac.dnase.peaks.sorted.bed


regions=results/chrombpnet/k562.merged.atac.dnase.peaks.bed
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir



