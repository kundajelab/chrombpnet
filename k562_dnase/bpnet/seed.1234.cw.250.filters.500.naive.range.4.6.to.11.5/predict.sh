#!/bin/bash

#fold to use for training 
fold=$1

#countweight
weight=$2

#gpu to use for training 
gpu=$3

#set seed for training
if [ -z "$4" ]
then
    seed=1234
else
    seed=$4
fi
echo "seed:$seed"

#output directory 
if [ -z "$5" ]
then
    outdir='.'
else
    outdir=$5
fi
echo "outdir:$outdir"


CUDA_VISIBLE_DEVICES=$gpu kerasAC_predict_tdb \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/dnase \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 2 \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 673 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --num_inputs 1 \
		    --num_outputs 2 \
		    --chrom_sizes ~/hg38.chrom.sizes \
		    --tiledb_stride 50 \
		    --fold $fold \
		    --genome hg38 \
		    --upsample_ratio_list_predict 1 \
		    --predictions_and_labels_hdf5 $outdir/predictions.k562.withdups.$seed\seed.$weight\counts.$fold \
		    --load_model_hdf5 $outdir/K562.profile.peaks.only.bpnet.withdups.$seed\seed.$weight\counts.$fold.hdf5 \
		    --tasks ENCSR000EOT \
		    --upsample_threads 1 \
		    --tdb_ambig_attribute ambig_peak \
		    --tdb_transformation_pseudocount 1
