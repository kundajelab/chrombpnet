#prefix=/srv/scratch/annashch/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/invivo_bias_model_step1
prefix=/data/chrombpnet/microglia/results/chrombpnet/ATAC/microglia/4_4_shifted_ATAC_10.11.2021_bias_filters_128/invivo_bias_model_step1
#CUDA_VISIBLE_DEVICES=3 kerasAC_form_modisco_inputs --model_hdf5 $prefix/microglia.0.hdf5 \
#		    --peak_file /srv/scratch/annashch/chrombpnet/microglia/data/microglia.idr.optimal.narrowPeak \
#		    --out_prefix $prefix/modisco/microglia.modisco \
#		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
#		    --chrom_sizes ~/hg38.chrom.sizes
kerasAC_run_modisco --hyp_imp $prefix/modisco/microglia.modisco.hyp.profile.npy \
		    --imp $prefix/modisco/microglia.modisco.observed.profile.npy  \
		    --onehot $prefix/modisco/microglia.modisco.seq.npy \
		    --outfile $prefix/modisco/microglia.modisco.output.profile \
		    --cores 40 \
		    --seqlets 25000
kerasAC_run_modisco --hyp_imp $prefix/modisco/microglia.modisco.hyp.count.npy \
		    --imp $prefix/modisco/microglia.modisco.observed.count.npy  \
		    --onehot $prefix/modisco/microglia.modisco.seq.npy \
		    --outfile $prefix/modisco/microglia.modisco.output.counts \
		    --cores 40 \
		    --seqlets 25000

