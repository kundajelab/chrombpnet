prefix=/oak/stanford/groups/akundaje/projects/encode_tier1_lines/SVM
cell=$1
fold=$2
for i in `seq 0 150`
do
    cat $prefix/$cell/genomewide_predictions/genomewidepredictions.$cell.$fold.$i >> $cell/svm_predictions_svmtrainset_genometestset/genomewidepredictions.$cell.$fold.all
    echo $i
done
