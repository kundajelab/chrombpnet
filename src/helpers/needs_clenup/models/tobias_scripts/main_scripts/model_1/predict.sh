#!/bin/bash
fold=$1

gpu=$2

#create a title for the model
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
tdb_array=$6
cell_line=$7
chrom_sizes=$8

echo "outdir:$outdir"
CUDA_VISIBLE_DEVICES=$gpu kerasAC_predict_tdb \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array  $tdb_array \
		    --tdb_partition_attribute_for_upsample idr_peak \
		    --tdb_partition_thresh_for_upsample 2 \
		    --tdb_partition_datasets_for_upsample $cell_line \
		    --tdb_input_source_attribute seq \
		    --tdb_input_aggregation None \
		    --tdb_input_transformation None \
		    --tdb_input_flank 1057 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation softmax None \
		    --num_inputs 1 \
		    --num_outputs 2 \
		    --chrom_sizes $chrom_sizes \
		    --tiledb_stride 50 \
		    --fold $fold \
		    --genome hg38 \
		    --upsample_ratio_list_predict 1 \
		    --predictions_and_labels_hdf5 $outdir/$model_name.$fold \
		    --json $outdir/$model_name.$fold.arch \
		    --weights $outdir/$model_name.$fold.weights \
		    --tdb_input_datasets seq \
		    --tdb_output_datasets $cell_line $cell_line \
		    --upsample_threads 1 \
		    --tdb_ambig_attribute ambig_peak \
