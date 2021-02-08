#!/bin/bash
task=$1
gpu=$2
prefix=/srv/scratch/annashch/5_cell_lines_bias_correction/svm
[ -d $task/dl_predictions_svmtrainset_genometestset ] || mkdir -p $task/dl_predictions_svmtrainset_genometestset
for fold in `seq 4 9`
do
    
    CUDA_VISIBLE_DEVICES=$gpu kerasAC_predict --index_data_path /srv/scratch/annashch/5_cell_lines_bias_correction/genomewide_labels/classificationlabels.SummitWithin200bpCenter.hdf5 \
			--input_data_path seq /srv/scratch/annashch/5_cell_lines_bias_correction/genomewide_gc/gc_hg38_nosmooth.hdf5 \
			--output_data_path /srv/scratch/annashch/5_cell_lines_bias_correction/genomewide_labels/classificationlabels.SummitWithin200bpCenter.hdf5 \
			--num_inputs 2 \
			--num_outputs 1 \
			--ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
			--load_model_hdf5 $prefix/$task/dl_model_svmtrainset/$task.$fold.classification.withgc.svmtrainset.hdf5 \
			--threads 20 \
			--max_queue_size 100 \
			--predictions_and_labels_hdf5 $prefix/$task/dl_predictions_svmtrainset_genometestset/$task.$fold.classification.withgc.dl.pred.svmtrainset.genometestset.hdf5 \
			--genome hg38 \
			--fold $fold \
			--batch_size 200 \
			--expand_dims \
			--tasks $task gc_fract \
			--index_task $task
done
