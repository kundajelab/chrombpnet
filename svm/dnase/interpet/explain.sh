#!/bin/bash
split=$1
task=$2
prefix=/oak/stanford/groups/akundaje/projects/enzymatic_bias_correction/svm/dnase/$task
model=$prefix/models/$task.$split.model.txt
gkmexplain  $prefix/fasta/fasta_interpret/$task.DNASE.svm.10kb.to.interpret.fa $model  $prefix/interpret/gkmexplain.$task.$split.dnase.txt
