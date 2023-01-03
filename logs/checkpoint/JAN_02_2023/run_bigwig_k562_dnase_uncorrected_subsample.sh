chrombpnet_nb=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
chrombpnet=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
bias=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/hint_atac.h5
celline=K562
gpu=3

output_dir=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/interpret/
mkdir $output_dir

#merge k562 peaks form atac and dnase
#atac_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/peaks.bed.gz
#dnase_peaks=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/caper/K562.ENFF205FNC/call-reproducibility_overlap/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/overlap.optimal_peak.narrowPeak.gz
#zcat $atac_peaks $dnase_peaks | bedtools sort -i stdin | uniq > results/chrombpnet/k562.merged.atac.dnase.peaks.sorted.bed


#regions=results/chrombpnet/k562.merged.atac.dnase.peaks.bed
regions=results/chrombpnet/DNASE_PE/K562/uncorrected_model_08.31.2022_filters_512_dil_8/uncorrected_model/interpret/K562.interpreted_regions_v2.bed 
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir



