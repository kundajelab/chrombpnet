fold=2
CUDA_VISIBLE_DEVICES=1 kerasAC_bpnet_shap_wrapper \
		    --model_hdf5 /srv/scratch/annashch/deeplearning/profile/k562_dnase/bpnet/seed.1234.cw.250.filters.500.naive.range.4.6.to.11.5/K562.profile.peaks.only.bpnet.withdups.1234seed.250counts.$fold.hdf5 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --bed_regions K562.dnase.idr.optimal_peak.narrowPeak.gz \
		    --bed_regions_center summit \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/dnase \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --tasks ENCSR000EOT \
		    --batch_size 200 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_input_source_attribute seq \
		    --tdb_input_flank 673 \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --out_pickle K562.DNASE.baseline.fold$fold.deepSHAP \
		    --num_threads 10
