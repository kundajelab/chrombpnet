#!/bin/bash
task=$1
fold=$2
prefix=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/svm/atac
#gkmtrain -m 10000 -v 2 -T 16 $prefix/$task/fasta/svm.inputs.$task.train.$fold.positives $prefix/$task/fasta/svm.inputs.$task.train.$fold.negatives $prefix/$task/models/model.$task.$fold 

gkmpredict -v 2 -T 16 $prefix/$task/fasta/svm.inputs.$task.test.$fold.positives $prefix/$task/models/model.$task.$fold.txt $prefix/$task/preds/$task.$fold.positives 
gkmpredict -v 2 -T 16 $prefix/$task/fasta/svm.inputs.$task.test.$fold.negatives $prefix/$task/models/model.$task.$fold.txt $prefix/$task/preds/$task.$fold.negatives 

