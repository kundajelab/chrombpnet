regions=results/chrombpnet/ATAC_PE/K562/data/peaks_no_blacklist.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/K562_wo_bias.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/K562_wo_bias_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig


regions=results/chrombpnet/ATAC_PE/K562/data/peaks_no_blacklist.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC/K562_unstranded.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/K562_unstranded_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

regions=results/chrombpnet/ATAC_PE/K562/data/peaks_no_blacklist.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/K562_w_bias.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/K562_w_bias_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.profile.bw
bigwig_stat=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.profile.stat
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.profile_scaled.bw

#python normalize_shap_bigwig_tracks.py --input_stat $bigwig_stat -bigwig $bigwig -o $obigwig

bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.counts.bw
bigwig_stat=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.counts.stat
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/ATAC_temp/merged.K562.counts_scaled.bw

#python normalize_shap_bigwig_tracks.py --input_stat $bigwig_stat -bigwig $bigwig -o $obigwig


#regions=results/chrombpnet/DNASE_PE/K562/data/peaks_no_blacklist.merged.bed
#bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_wo_bias.bw
#obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_wo_bias_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

regions=results/chrombpnet/DNASE_PE/K562/data/peaks_no_blacklist.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_w_bias.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_w_bias_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig


bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.profile.bw
bigwig_stat=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.profile.stat
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.profile_scaled.bw

#python normalize_shap_bigwig_tracks.py --input_stat $bigwig_stat -bigwig $bigwig -o $obigwig

bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.counts.bw
bigwig_stat=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.counts.stat
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/merged.K562.counts_scaled.bw

#python normalize_shap_bigwig_tracks.py --input_stat $bigwig_stat -bigwig $bigwig -o $obigwig

bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_unstranded.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/K562/DNASE/K562_unstranded_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig

