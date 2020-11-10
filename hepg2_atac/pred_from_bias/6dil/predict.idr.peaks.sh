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
echo "outdir:$outdir"
CUDA_VISIBLE_DEVICES=$gpu kerasAC_predict_tdb \
		    --bed_regions HEPG2.atac.idr.optimal_peak.narrowPeak.gz \
		    --center_on_summit \
		    --batch_size 25 \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --tdb_array /srv/scratch/annashch/encode_dnase_tiledb/db/atac \
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
		    --genome hg38 \
		    --predictions_and_labels_hdf5 $outdir/$model_name.$fold \
		    --load_model_hdf5 $model_name.$fold.hdf5 \
		    --tasks HEPG2 \
		    --tdb_transformation_pseudocount 0.001
