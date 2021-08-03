#!/bin/bash

#fold to use for training 
fold=$1

#gpu to use for training 
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
params=$6
tdb_array=$7
cell_line=$8
arch_file=$9

echo "outdir:$outdir"
CUDA_VISIBLE_DEVICES=$gpu kerasAC_train \
		    --seed $seed \
		    --batch_size 20 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array $tdb_array \
		    --tdb_partition_attribute_for_upsample overlap_peak \
		    --tdb_partition_thresh_for_upsample 1 \
		    --tdb_partition_datasets_for_upsample $cell_line \
		    --tdb_input_source_attribute control_count_bigwig_unstranded_5p control_count_bigwig_unstranded_5p \
		    --tdb_input_aggregation None sum \
		    --tdb_input_transformation None log \
		    --tdb_input_flank 500 500 \
		    --tdb_output_source_attribute count_bigwig_unstranded_5p count_bigwig_unstranded_5p \
		    --tdb_output_flank 500 500 \
		    --tdb_output_aggregation None sum \
		    --tdb_output_transformation None log \
		    --tdb_ambig_attribute ambig_peak \
		    --tdb_input_min None None \
		    --tdb_input_max None None \
		    --tdb_output_min None 4.0 \
		    --tdb_output_max None 11.5 \
		    --num_inputs 2 \
		    --num_outputs 2 \
		    --fold $fold \
		    --genome hg38 \
		    --num_train 100000 \
		    --num_valid 10000 \
		    --num_tasks 1 \
		    --upsample_threads 24 \
		    --threads 0 \
		    --max_queue_size 20 \
		    --patience 3 \
		    --patience_lr 2 \
		    --model_prefix $outdir/$model_name.$fold \
		    --architecture_from_file $arch_file \
		    --model_params $params \
		    --tdb_input_datasets $cell_line $cell_line \
		    --tdb_output_datasets $cell_line $cell_line \
		    --upsample_ratio_list_train 1.0 \
		    --upsample_ratio_list_eval 1.0 \
		    --trackables logcount_predictions_loss loss profile_predictions_loss val_logcount_predictions_loss val_loss val_profile_predictions_loss

