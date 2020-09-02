#!/bin/bash
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1

#trained on peaks
#CUDA_VISIBLE_DEVICES=2 kerasAC_predict \
#		    --batch_size 20 \
#		    --ref_fasta /mnt/data/annotations/by_release/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
#		    --tdb_indexer tasks.tsv \
#		    --tdb_inputs seq \
#		    --tdb_input_source_attribute seq \
#		    --tdb_input_aggregation None \
#		    --tdb_input_transformation None \
#		    --tdb_input_flank 6500 \
#		    --tdb_outputs tasks.tsv \
#		    --tdb_output_source_attribute fc_bigwig \
#		    --tdb_output_flank 1585 \
#		    --tdb_output_aggregation None \
#		    --tdb_output_transformation asinh \
#		    --num_inputs 1 \
#		    --num_outputs 1 \
#		    --chrom_sizes /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
#		    --predict_chroms chr1 \
#		    --threads 10 \
#		    --max_queue_size 500 \
#		    --model_hdf5 ATAC.pseudobulk.ADPD.Cluster24.profile.peaks.only.0 \
#		    --predictions_hdf5 predictions.genomewide.microglia.0.chr1.trainedonpeaks \
#		    --tiledb_stride 10000

#trained genomewide
CUDA_VISIBLE_DEVICES=2 kerasAC_predict \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/annotations/by_release/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_indexer tasks.tsv \
		    --tdb_inputs seq \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 6500 \
		    --tdb_outputs tasks.tsv \
		    --tdb_output_source_attribute fc_bigwig \
		    --tdb_output_flank 1585 \
		    --tdb_output_aggregation None \
		    --tdb_output_transformation asinh \
		    --num_inputs 1 \
		    --num_outputs 1 \
		    --chrom_sizes /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
		    --predict_chroms chr1 \
		    --threads 10 \
		    --max_queue_size 500 \
		    --model_hdf5 ATAC.pseudobulk.ADPD.Cluster24.profile.0 \
		    --predictions_hdf5 predictions.genomewide.microglia.0.chr1 \
		    --tiledb_stride 10000
