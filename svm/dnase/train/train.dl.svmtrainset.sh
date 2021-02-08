#!/bin/bash
task=$1
gpu=$2
prefix=/srv/scratch/annashch/5_cell_lines_bias_correction/svm/
[ -d $task/dl_model_svmtrainset ] || mkdir -p $task/dl_model_svmtrainset
for fold in `seq 0 9`
do    
    CUDA_VISIBLE_DEVICES=$gpu kerasAC_train --index_data_path $prefix/$task/$task.dl.inputs.train.$fold.$task.labels.hdf5 \
		    --input_data_path seq $prefix/$task/$task.dl.inputs.train.$fold.$task.gc.hdf5 \
		    --output_data_path  $prefix/$task/$task.dl.inputs.train.$fold.$task.labels.hdf5 \
		    --num_inputs 2 \
		    --num_outputs 1 \
		    --model_prefix $prefix/$task/dl_model_svmtrainset/$task.$fold.classification.withgc.svmtrainset \
		    --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
		    --batch_size 200 \
		    --architecture_spec functional_basset_classification_gc_corrected \
		    --num_train 100000 \
		    --num_valid 10000 \
		    --num_tasks 1 \
		    --threads 20 \
		    --max_queue_size 100 \
		    --init_weights /srv/scratch/annashch/deeplearning/encode-roadmap.dnase_tf-chip.batch_256.params.npz \
		    --patience 3 \
		    --patience_lr 2 \
		    --expand_dims \
		    --tasks $task gc \
		    --index_tasks $task \
		    --fold $fold \
		    --genome hg38
done
