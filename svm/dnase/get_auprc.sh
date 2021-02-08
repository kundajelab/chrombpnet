for task in GM12878 H1ESC HEPG2 IMR90 K562
do
    for fold in `seq 0 9`
    do
	python get_auprc.py $task $fold  > /srv/scratch/annashch/5_cell_lines_bias_correction/svm/perf.svm.svmtrainset.genometestset/perf.txt.$task.$fold & 
    done
done
