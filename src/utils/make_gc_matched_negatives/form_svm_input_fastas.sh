task=$1
ref=$2
python $PWD/main_scripts/make_gc_matched_negatives/form_svm_input_fastas.py --outf $task/bpnet.inputs.test.0 $task/bpnet.inputs.train.0 \
       --neg_pickle $task/candidate.negatives.gc.p \
       --overwrite_outf \
       --ref_fasta $ref \
       --peaks $task/bpnet.peaks.test.0.gc.seq $task/bpnet.peaks.train.0.gc.seq 
