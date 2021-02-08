
for task in IMR90  #GM12878 H1ESC HEPG2 IMR90 K562
do
    #DELETE THE OLD ONES!
    rm $task/$task.svm.inputs.t* 
    for chrom in `seq 1 22` X Y
    do
	python form_svm_input_fastas.py --outf $task/$task.svm.inputs.test.0  $task/$task.svm.inputs.test.4  $task/$task.svm.inputs.test.8   $task/$task.svm.inputs.train.2  $task/$task.svm.inputs.train.6 $task/$task.svm.inputs.test.1  $task/$task.svm.inputs.test.5  $task/$task.svm.inputs.test.9   $task/$task.svm.inputs.train.3  $task/$task.svm.inputs.train.7 $task/$task.svm.inputs.test.2  $task/$task.svm.inputs.test.6  $task/$task.svm.inputs.train.0  $task/$task.svm.inputs.train.4  $task/$task.svm.inputs.train.8 $task/$task.svm.inputs.test.3  $task/$task.svm.inputs.test.7  $task/$task.svm.inputs.train.1  $task/$task.svm.inputs.train.5  $task/$task.svm.inputs.train.9 \
	       --chrom chr$chrom \
	       --neg_pickle $task/$task.chr$chrom.negative.gc.seq.pickle \
	       --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	       --peaks $task/$task.svm.peaks.test.0.gc.seq  $task/$task.svm.peaks.test.4.gc.seq  $task/$task.svm.peaks.test.8.gc.seq   $task/$task.svm.peaks.train.2.gc.seq  $task/$task.svm.peaks.train.6.gc.seq $task/$task.svm.peaks.test.1.gc.seq  $task/$task.svm.peaks.test.5.gc.seq  $task/$task.svm.peaks.test.9.gc.seq   $task/$task.svm.peaks.train.3.gc.seq  $task/$task.svm.peaks.train.7.gc.seq $task/$task.svm.peaks.test.2.gc.seq  $task/$task.svm.peaks.test.6.gc.seq  $task/$task.svm.peaks.train.0.gc.seq  $task/$task.svm.peaks.train.4.gc.seq  $task/$task.svm.peaks.train.8.gc.seq $task/$task.svm.peaks.test.3.gc.seq  $task/$task.svm.peaks.test.7.gc.seq  $task/$task.svm.peaks.train.1.gc.seq  $task/$task.svm.peaks.train.5.gc.seq  $task/$task.svm.peaks.train.9.gc.seq &
    done
done


