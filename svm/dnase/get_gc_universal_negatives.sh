for task in GM12878 HEPG2 H1ESC IMR90 K562
do
    for i in `seq 1 22` X Y
    do
	python /srv/scratch/annashch/5_cell_lines_bias_correction/genomewide_gc/get_gc_content.py \
	       --input_bed $task/$task.chr$i.universal_negatives.txt \
	       --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	       --out_prefix $task/$task.chr$i.negative.gc.seq \
	       --store_seq &
    done
done
