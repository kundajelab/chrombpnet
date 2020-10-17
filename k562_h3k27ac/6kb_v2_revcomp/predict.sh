fold=0
out_dir='.'
model_prefix=tmp
CUDA_VISIBLE_DEVICES=0 kerasAC_predict_tdb \
		    --batch_size 100 \
		    --ref_fasta /users/annashch/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/histone \
		    --num_inputs 3 \
		    --num_outputs 2 \
		    --tdb_input_datasets seq K562,K562 K562,K562 \
		    --tdb_output_datasets K562,K562 K562,K562 \
		    --tdb_input_source_attribute seq control_count_bigwig_plus_5p,control_count_bigwig_minus_5p control_count_bigwig_plus_5p,control_count_bigwig_minus_5p \
		    --tdb_output_source_attribute count_bigwig_plus_5p,count_bigwig_minus_5p count_bigwig_plus_5p,count_bigwig_minus_5p \
		    --tdb_input_flank 3000 500,500 500,500 \
		    --tdb_output_flank 500,500 500,500 \
		    --tdb_input_min None None,None None,None \
		    --tdb_output_min None,None None,None \
		    --tdb_input_max None None,None None,None \
		    --tdb_output_max None,None None,None \
		    --tdb_input_aggregation None None,None sum,sum \
		    --tdb_input_transformation None None,None log,log \
		    --tdb_output_aggregation None,None sum,sum \
		    --tdb_output_transformation None,None log,log \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --load_model_hdf5 tmp.hdf5 \
		    --bed_regions overlap.optimal_peak.narrowPeak.gz \
		    --bed_regions_center summit \
		    --predictions_and_labels_hdf5 preds.hdf5
