for cell in GM12878 H1ESC HEPG2 IMR90 K562
do
    for fold in `seq 0 9`
    do
	rm $cell/svm_predictions_svmtrainset_genometestset/genomewidepredictions.$cell.$fold.all
	./collapse.sh $cell $fold & 
    done
done
