#!/bin/bash

#fold to use for training 
fold=$1


#gpu to use for training 
gpu=$2

#get model name
model_name=$3

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
		    --batch_size 100 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array /srv/scratch/annashch/deeplearning/profile/gata/db/gata2 \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 2 \
		    --tdb_input_source_attribute seq control_count_bigwig_plus_5p,control_count_bigwig_minus_5p control_count_bigwig_plus_5p,control_count_bigwig_minus_5p \
		    --tdb_input_aggregation None None sum \
		    --tdb_input_transformation None None log \
		    --tdb_input_flank 673 500 500 \
		    --tdb_output_source_attribute count_bigwig_plus_5p,count_bigwig_minus_5p count_bigwig_plus_5p,count_bigwig_minus_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
     		    --num_inputs 3 \
		    --num_outputs 2 \
		    --tdb_ambig_attribute ambig_peak \
		    --chrom_sizes ~/hg38.chrom.sizes \
		    --fold $fold \
		    --genome hg38 \
		    --upsample_ratio_list_predict 1 \
		    --predictions_and_labels_hdf5 $outdir/$model_name.$fold \
		    --load_model_hdf5 $outdir/$model_name.$fold.hdf5 \
		    --tasks GATA2 \
		    --upsample_threads 1 \
		    --tdb_transformation_pseudocount 0.001
