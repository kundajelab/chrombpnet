#!/bin/bash
#--train_chroms chr2 chr3 chr4 chr5 chr6 chr7 chr9 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
#--validation_chroms chr8 chr10 \
#--train_chroms chr21 \
#--validation_chroms chr22 \

#export OMP_NUM_THREADS=1
#export USE_SIMPLE_THREADED_LEVEL3=1
fold=0
CUDA_VISIBLE_DEVICES=0 kerasAC_predict_tdb \
		    --batch_size 20 \
		    --ref_fasta GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_indexer tasks.nodups.tsv \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 2 \
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
		    --chrom_sizes hg38.chrom.sizes \
		    --tiledb_stride 50 \
		    --fold $fold \
		    --genome hg38 \
		    --upsample_ratio_list_predict 1 \
		    --predictions_and_labels_hdf5 predictions.k562.nodups.0 \
		    --load_model_hdf5 K562.profile.peaks.only.bpnet.nodups.0.hdf5
