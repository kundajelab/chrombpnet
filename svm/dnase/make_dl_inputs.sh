#generate training splits 
for task in GM12878 H1ESC HEPG2 IMR90 K562
do
    for split in `seq 0 9`
    do
	for purpose in train test
		       
	do
	    python make_dl_inputs.py --tasks $task \
	       --positives $task/$task.svm.inputs.$purpose.$split.positives.bed \
	       --negatives $task/$task.svm.inputs.$purpose.$split.negatives.bed \
	       --store_gc \
	       --out_prefix $task/$task.dl.inputs.$purpose.$split
	done
    done
done
