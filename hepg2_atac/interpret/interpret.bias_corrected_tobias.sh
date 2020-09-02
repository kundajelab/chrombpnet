fold=2
CUDA_VISIBLE_DEVICES=0 kerasAC_bpnet_shap_wrapper \
		    --model_hdf5 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_tobias/hepg2.atac.with.tobias.bias.$fold.hdf5 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --bed_regions HEPG2.atac.idr.optimal_peak.narrowPeak.gz \
		    --bed_regions_center summit \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/atac \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --tasks HEPG2 \
		    --batch_size 200 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_input_source_attribute seq \
		    --tdb_input_flank 673 \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --out_pickle HEPG2.ATAC.bias_corrected_tobias.fold$fold.deepSHAP \
		    --num_threads 10
fold=3
CUDA_VISIBLE_DEVICES=0 kerasAC_bpnet_shap_wrapper \
		    --model_hdf5 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_tobias/hepg2.atac.with.tobias.bias.$fold.hdf5 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --bed_regions HEPG2.atac.idr.optimal_peak.narrowPeak.gz \
		    --bed_regions_center summit \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/atac \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --tasks HEPG2 \
		    --batch_size 200 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_input_source_attribute seq \
		    --tdb_input_flank 673 \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --out_pickle HEPG2.ATAC.bias_corrected_tobias.fold$fold.deepSHAP \
		    --num_threads 10

fold=4
CUDA_VISIBLE_DEVICES=0 kerasAC_bpnet_shap_wrapper \
		    --model_hdf5 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/with_bias_tobias/hepg2.atac.with.tobias.bias.$fold.hdf5 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --bed_regions HEPG2.atac.idr.optimal_peak.narrowPeak.gz \
		    --bed_regions_center summit \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/atac \
		    --chrom_sizes /users/annashch/hg38.chrom.sizes \
		    --tasks HEPG2 \
		    --batch_size 200 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_input_source_attribute seq \
		    --tdb_input_flank 673 \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --out_pickle HEPG2.ATAC.bias_corrected_tobias.fold$fold.deepSHAP \
		    --num_threads 10
