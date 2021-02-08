#!/bin/bash
split=$1
task=$2
prefix=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/svm/atac/$task
model=$prefix/models/model.$task.$split.txt
gkmexplain  $prefix/fasta/fasta_interpret/$task.ATAC.svm.10kb.to.interpret.fa $model  $prefix/interpret/gkmexplain.$task.$split.atac.txt
