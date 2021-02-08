#for task in GM12878 HEPG2 H1ESC IMR90 K562
for task in HEPG2 H1ESC IMR90
do
    python get_svm_peak_splits.py \
	   --narrowPeak $task/$task.idr.narrowPeak.gz \
	   --ntrain 60000 \
	   --out_prefix $task/$task.svm.peaks \
	   --genome hg38
done
