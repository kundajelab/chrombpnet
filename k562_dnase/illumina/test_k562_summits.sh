#!/bin/bash
#--train_chroms chr2 chr3 chr4 chr5 chr6 chr7 chr9 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
#--validation_chroms chr8 chr10 \
#--train_chroms chr21 \
#--validation_chroms chr22 \

#export OMP_NUM_THREADS=1
#export USE_SIMPLE_THREADED_LEVEL3=1
CUDA_VISIBLE_DEVICES=3 kerasAC_predict \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_indexer tasks.tsv \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 2 \
		    --tdb_inputs seq \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 6500 \
		    --tdb_outputs tasks.tsv \
		    --tdb_output_source_attribute fc_bigwig \
		    --tdb_output_flank 1500 \
		    --tdb_output_aggregation None \
		    --num_inputs 1 \
		    --num_outputs 1 \
		    --tdb_output_transformation asinh \
		    --chrom_sizes /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
		    --tiledb_stride 50 \
		    --predict_chroms chr1 \
		    --upsample_ratio_list_predict 1 \
		    --predictions_and_labels_hdf5 predictions.k562.0.chr1 \
		    --threads 10 \
		    --max_queue_size 10 \
		    --model_hdf5 K562.profile.peaks.only.illumina.0
