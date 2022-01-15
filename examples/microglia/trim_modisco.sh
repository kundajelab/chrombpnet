### BIAS MODELS ### 
bias_modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/invivo_bias_model_step1/modisco
mkdir $bias_modisco_folder/microglia_profile_trimmed_modisco_plots
mkdir $bias_modisco_folder/microglia_counts_trimmed_modisco_plots

kerasAC_trim_modisco --modisco_hits $bias_modisco_folder/microglia.modisco.output.profile --trim_thresh 0.45 --trim_extend 1
mv *png $bias_modisco_folder/microglia_profile_trimmed_modisco_plots
kerasAC_trim_modisco --modisco_hits $bias_modisco_folder/microglia.modisco.output.counts --trim_thresh 0.45 --trim_extend 1
mv *png $bias_modisco_folder/microglia_counts_trimmed_modisco_plots

### MODISCO MODELS ###
modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/final_model_step3/unplug/modisco
mkdir $modisco_folder/microglia_profile_trimmed_modisco_plots
mkdir $modisco_folder/microglia_counts_trimmed_modisco_plots

kerasAC_trim_modisco --modisco_hits $modisco_folder/microglia.modisco.output.profile --trim_thresh 0.45 --trim_extend 1
mv *png $modisco_folder/microglia_profile_trimmed_modisco_plots
kerasAC_trim_modisco --modisco_hits $modisco_folder/microglia.modisco.output.counts --trim_thresh 0.45 --trim_extend 1
mv *png $modisco_folder/microglia_counts_trimmed_modisco_plots
