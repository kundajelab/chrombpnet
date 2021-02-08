#!/bin/bash
task=$1
fold=$2
gkmtrain -m 10000 -v 2 -T 16 $task/$task.svm.inputs.train.$fold.positives $task/$task.svm.inputs.train.$fold.negatives $task/model.$task.$fold
