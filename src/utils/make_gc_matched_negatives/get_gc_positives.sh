task=$1
ref_fasta=$2
flank_size=$3
for split in 0
do
    for dataset in train test
    do
	python $PWD/src/utils/make_gc_matched_negatives/get_gc_content.py \
	       --input_bed $task/bpnet.peaks.$dataset.$split \
	       --ref_fasta $ref_fasta \
	       --out_prefix $task/bpnet.peaks.$dataset.$split.gc.seq \
	       --center_summit \
	       --flank_size $flank_size \
	       --store_seq
    done
done
