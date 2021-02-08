for task in GM12878 H1ESC HEPG2 IMR90 K562
do
    for fold in `seq 0 9`
    do
	python pad_na.py $task $fold & 
    done
done
