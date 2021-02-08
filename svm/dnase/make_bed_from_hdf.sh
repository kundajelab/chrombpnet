#create bed files from hdf5 dl model genomewide predictions
for task in GM12878 H1ESC HEPG2 IMR90 K562
do
    for fold in `seq 0 9`
    do
	python make_bed_from_hdf.py $task $fold &
    done
done
