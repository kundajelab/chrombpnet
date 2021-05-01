for fold in 0 
do
    CUDA_VISIBLE_DEVICES=0 kerasAC_bpnet_shap_wrapper \
			--model_hdf5 imr90.dnase.with.bias.unplugged.0.hdf5 \
			--ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
			--bed_regions IMR90.dnase.idr.optimal_peak.narrowPeak.gz \
			--bed_regions_center summit \
			--tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/dnase \
			--chrom_sizes /users/annashch/hg38.chrom.sizes \
			--tdb_input_datasets seq \
			--tdb_output_datasets ENCSR477RTP ENCSR477RTP \
			--batch_size 200 \
			--tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
			--tdb_output_flank 500 500 \
			--tdb_output_aggregation None sum \
			--tdb_output_transformation None log \
			--tdb_input_source_attribute seq \
			--tdb_input_flank 1057 \
			--tdb_input_aggregation None \
			--tdb_input_transformation None \
			--out_pickle IMR90.DNASE.bias_corrected_bpnet_tobias.fold$fold.deepSHAP \
			--num_threads 20 \
			--num_inputs 1 \
			--num_outputs 2 \
			--task_index 0
done
