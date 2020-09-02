#!/bin/bash
#export OMP_NUM_THREADS=1
#export USE_SIMPLE_THREADED_LEVEL3=1
fold=0
CUDA_VISIBLE_DEVICES=0 kerasAC_train \
		    --batch_size 20 \
		    --ref_fasta GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_indexer tasks.nodups.tsv \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 1 \
		    --tdb_inputs seq \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 673 \
		    --tdb_outputs tasks.nodups.tsv tasks.nodups.tsv \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None asinh \
		    --num_inputs 1 \
		    --num_outputs 2 \
		    --fold $fold \
		    --genome hg38 \
		    --chrom_sizes /mnt/data/annotations/by_release/hg38/hg38.chrom.sizes \
		    --upsample_ratio_list_train 1 \
		    --upsample_ratio_list_eval 1 \
		    --num_train 100000 \
		    --num_valid 10000 \
		    --num_tasks 1 \
		    --threads 0 \
		    --max_queue_size 20 \
		    --patience 3 \
		    --patience_lr 2 \
		    --model_prefix K562.profile.peaks.only.bpnet.nodups.$fold \
		    --architecture_spec profile_bpnet_dnase \
		    --use_multiprocessing False
