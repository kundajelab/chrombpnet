for task in GM12878 HEPG2 H1ESC IMR90 K562
do
    for chrom in `seq 1 22` X Y
    do
	python get_negatives.py $task chr$chrom &
    done
done
