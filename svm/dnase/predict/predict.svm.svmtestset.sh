#!/bin/bash
task=$1
fold=$2
gkmtrain -m 10000 -v 2 -T 16 $task/$task.svm.inputs.train.$fold.positives $task/$task.svm.inputs.train.$fold.negatives $task/model.$task.$fold
gkmpredict -v 2 -T 16 $task/$task.svm.inputs.test.$fold.positives $task/model.$task.$fold $task/predictions.$task.$fold.positives
gkmpredict -v 2 -T 16 $task/$task.svm.inputs.test.$fold.negatives $task/model.$task.$fold $task/predictions.$task.$fold.negatives    
