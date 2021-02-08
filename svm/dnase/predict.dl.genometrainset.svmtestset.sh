#!/bin/bash
task=$1
gpu=$2
prefix=/srv/scratch/annashch/5_cell_lines_bias_correction/svm
[ -d $task/dl_predictions_genometrainset_svmtestset ] || mkdir -p $task/dl_predictions_genometrainset_svmtestset
for fold in `seq 0 9`
do
    
    CUDA_VISIBLE_DEVICES=$gpu kerasAC_predict --index_data_path $prefix/$task/$task.dl.inputs.test.$fold.$task.labels.hdf5 \
			--input_data_path seq $prefix/$task/$task.dl.inputs.test.$fold.$task.gc.hdf5 \
			--output_data_path $prefix/$task/$task.dl.inputs.test.$fold.$task.labels.hdf5 \
			--num_inputs 2 \
			--num_outputs 1 \
			--ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
			--load_model_hdf5 /srv/scratch/annashch/5_cell_lines_bias_correction/gc_covariate/classification/$task/DNASE.$task.classificationlabels.withgc.$fold \
			--threads 20 \
			--max_queue_size 100 \
			--predictions_and_labels_hdf5  $prefix/$task/dl_predictions_genometrainset_svmtestset/$task.$fold.classification.withgc.dl.pred.genometrainset.svmtestset.hdf5 \
			--genome hg38 \
			--fold $fold \
			--batch_size 1000 \
			--expand_dims \
			--tasks $task gc
done
