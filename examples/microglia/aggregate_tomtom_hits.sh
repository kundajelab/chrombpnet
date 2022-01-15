bias_modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/invivo_bias_model_step1/modisco
modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/final_model_step3/unplug/modisco

kerasAC_aggregate_tomtom_motifs --tomtom_results $bias_modisco_folder/tomtom.bias.microglia.modisco.output.profile \
				--modisco_meme_file_dir $bias_modisco_folder/profile_patterns \
				--outf $bias_modisco_folder/microglia.bias.sig.meme.aggregate.profile.txt
kerasAC_aggregate_tomtom_motifs --tomtom_results $bias_modisco_folder/tomtom.bias.microglia.modisco.output.counts \
				--modisco_meme_file_dir $bias_modisco_folder/count_patterns \
				--outf $bias_modisco_folder/microglia.bias.sig.meme.aggregate.counts.txt


kerasAC_aggregate_tomtom_motifs --tomtom_results $modisco_folder/tomtom.microglia.modisco.output.profile \
				--modisco_meme_file_dir $modisco_folder/profile_patterns \
				--outf $modisco_folder/microglia.sig.meme.aggregate.profile.txt
kerasAC_aggregate_tomtom_motifs --tomtom_results $modisco_folder/tomtom.microglia.modisco.output.counts \
				--modisco_meme_file_dir $modisco_folder/count_patterns \
				--outf $modisco_folder/microglia.sig.meme.aggregate.counts.txt

