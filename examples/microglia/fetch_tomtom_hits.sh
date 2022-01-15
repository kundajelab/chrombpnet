bias_modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/invivo_bias_model_step1/modisco
modisco_folder=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/final_model_step3/unplug/modisco

mkdir $bias_modisco_folder/profile_patterns
mkdir $bias_modisco_folder/count_patterns
mkdir $modisco_folder/profile_patterns
mkdir $modisco_folder/count_patterns

kerasAC_tomtom_scan_modisco -m $bias_modisco_folder/microglia.modisco.output.profile \
			    -o $bias_modisco_folder/tomtom.bias.microglia.modisco.output.profile \
			    -d motifs.meme.txt \
			    -n 5 \
			    -th 0.45 \
			    -q 0.1
mv microglia.*.txt $bias_modisco_folder/profile_patterns/


kerasAC_tomtom_scan_modisco -m $bias_modisco_folder/microglia.modisco.output.counts \
			    -o $bias_modisco_folder/tomtom.bias.microglia.modisco.output.counts \
			    -d motifs.meme.txt \
			    -n 5 \
			    -th 0.45 \
			    -q 0.1
mv microglia.*.txt $bias_modisco_folder/count_patterns/



kerasAC_tomtom_scan_modisco -m $modisco_folder/microglia.modisco.output.profile \
			    -o $modisco_folder/tomtom.microglia.modisco.output.profile \
			    -d motifs.meme.txt \
			    -n 5 \
			    -th 0.45 \
			    -q 0.1
mv microglia*txt $modisco_folder/profile_patterns/


kerasAC_tomtom_scan_modisco -m $modisco_folder/microglia.modisco.output.counts \
			    -o $modisco_folder/tomtom.microglia.modisco.output.counts \
			    -d motifs.meme.txt \
			    -n 5 \
			    -th 0.45 \
			    -q 0.1 
mv microglia*txt $modisco_folder/count_patterns/

