for task in GM12878 HEPG2 H1ESC IMR90 K562
do
    for i in `seq 1 22` X Y
    do
	python pickle_negatives.py $task chr$i &
    done
done
