mkdir results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret/tfhits/
python src/evaluation/invivo_footprints/tf_modiscohits.py --outdir=results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret/tfhits/ \
	results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret/K562.profile_scores.h5 \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/K562/ATAC_PE_12.30.2021/SIGNAL/modisco/modisco_results_allChroms_profile.hdf5 \
	results/chrombpnet/ATAC_PE/K562/ATAC_PE_12.30.2021/chrombpnet_model/interpret/K562.interpreted_regions.bed \
