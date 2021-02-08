for task in GM12878 H1ESC HEPG2 IMR90 K562
do
    for f in `ls $task/$task.svm.inputs.t*`
    do
	echo $f
	grep ">" $f | sed --expression='s/>//g' | sed --expression='s/\_/\t/g' > $f.bed
    done
done
