#!/bin/bash
#export OMP_NUM_THREADS=1
#export USE_SIMPLE_THREADED_LEVEL3=1

CUDA_VISIBLE_DEVICES=0 kerasAC_train \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_indexer tasks.tsv \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 1 \
		    --tdb_inputs seq \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 6500 \
		    --tdb_outputs tasks.tsv \
		    --tdb_output_source_attribute fc_bigwig \
		    --tdb_output_flank 1500 \
		    --tdb_output_aggregation None \
		    --tdb_output_transformation asinh \
		    --num_inputs 1 \
		    --num_outputs 1 \
		    --train_chroms chr2 chr3 chr4 chr5 chr6 chr7 chr9 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
		    --validation_chroms chr8 chr10 \
		    --chrom_sizes /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
		    --upsample_ratio_list_train 1 \
		    --upsample_ratio_list_eval 0 \
		    --num_train 100000 \
		    --num_valid 10000 \
		    --num_tasks 1 \
		    --threads 0 \
		    --max_queue_size 20 \
		    --patience 3 \
		    --patience_lr 2 \
		    --model_hdf5 K562.profile.peaks.only.illumina.0 \
		    --architecture_spec profile_illumina \
		    --use_multiprocessing False
