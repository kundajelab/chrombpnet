fold=0
CUDA_VISIBLE_DEVICES=3 python profile_deepshap_for_shuffled_seq.py \
		    --model_hdf5 /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.$fold.hdf5 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --bed_regions K562.dnase.idr.optimal_peak.narrowPeak \
		    --bed_regions_center summit \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/dnase \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --tdb_input_datasets seq \
		    --tdb_output_datasets ENCSR000EOT ENCSR000EOT \
		    --batch_size 200 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_input_source_attribute seq \
		    --tdb_input_flank 673 \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --out_pickle K562.DNASE.bias_corrected_bpnet_tobias.unplugged.SHUFFLED.fold$fold.deepSHAP \
		    --num_threads 20 \
		    --num_inputs 1 \
		    --num_outputs 2

