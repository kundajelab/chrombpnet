task=$1
idr=$2
python $PWD/main_scripts/make_gc_matched_negatives/get_svm_peak_splits.py \
       --narrowPeak $idr \
       --folds 0 \
       --out_prefix $task/bpnet.peaks \
       --genome hg38
