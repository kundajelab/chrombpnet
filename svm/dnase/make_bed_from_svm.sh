for cell in GM12878 #H1ESC HEPG2 IMR90 K562
do
    for fold in 6 #`seq 0 9`
    do
	
	python make_bed_from_svm.py $cell/svm_predictions_svmtrainset_genometestset/genomewidepredictions.$cell.$fold.all 
    done
done
