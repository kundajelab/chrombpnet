#for task in GM12878 HEPG2 IMR90 H1ESC
#do
#    idr=$task.atac.idr.optimal_peak.narrowPeak.gz
#
#    python ~/SVM_pipelines/make_inputs/get_svm_peak_splits.py \
#	   --narrowPeak $task/$idr \
#	   --ntrain 60000 \
#	   --out_prefix $task/svm.peaks.$task \
#	   --genome hg38 &
#done

#for task in GM12878 HEPG2 IMR90 H1ESC
#do
#    for split in `seq 0 9`
#    do
#	for dataset in train test
#	do
#	    python ~/SVM_pipelines/make_inputs/get_gc_content.py \
#		   --input_bed $task/svm.peaks.$task.$dataset.$split \
#		   --ref_fasta /mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
#		   --out_prefix $task/svm.peaks.$task.$dataset.$split.gc.seq \
#		   --center_summit \
#		   --flank_size 500 \
#		   --store_seq  &
#	done
#    done
#done

#get the inverse intersection of idr peak file and all gc genome bins
#for task in GM12878 HEPG2 IMR90 H1ESC
#do
#    idr=$task.atac.idr.optimal_peak.narrowPeak.gz
#    bedtools intersect -v -a /srv/scratch/annashch/bias_correction/genomewide_gc/gc_hg38_nosmooth.tsv -b $task/$idr > $task/$task.candidate.negatives.tsv &
#done


for task in GM12878 HEPG2 IMR90 H1ESC
do
    python ~/SVM_pipelines/make_inputs/get_chrom_gc_region_dict.py --input_bed $task/$task.candidate.negatives.tsv --outf $task/$task.candidate.negatives.gc.p 
done

for task in GM12878 HEPG2 IMR90 H1ESC
do
    ref=/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    python ~/SVM_pipelines/make_inputs/form_svm_input_fastas.py --outf $task/svm.inputs.$task.test.0 $task/svm.inputs.$task.test.1 $task/svm.inputs.$task.test.2 $task/svm.inputs.$task.test.3 $task/svm.inputs.$task.test.4 $task/svm.inputs.$task.test.5 $task/svm.inputs.$task.test.6 $task/svm.inputs.$task.test.7 $task/svm.inputs.$task.test.8 $task/svm.inputs.$task.test.9 $task/svm.inputs.$task.train.0 $task/svm.inputs.$task.train.1 $task/svm.inputs.$task.train.2 $task/svm.inputs.$task.train.3 $task/svm.inputs.$task.train.4 $task/svm.inputs.$task.train.5 $task/svm.inputs.$task.train.6 $task/svm.inputs.$task.train.7 $task/svm.inputs.$task.train.8 $task/svm.inputs.$task.train.9 \
	   --neg_pickle $task/$task.candidate.negatives.gc.p \
	   --overwrite_outf \
	   --ref_fasta $ref \
	   --peaks $task/svm.peaks.$task.test.0.gc.seq $task/svm.peaks.$task.test.1.gc.seq $task/svm.peaks.$task.test.2.gc.seq $task/svm.peaks.$task.test.3.gc.seq $task/svm.peaks.$task.test.4.gc.seq $task/svm.peaks.$task.test.5.gc.seq $task/svm.peaks.$task.test.6.gc.seq $task/svm.peaks.$task.test.7.gc.seq $task/svm.peaks.$task.test.8.gc.seq $task/svm.peaks.$task.test.9.gc.seq $task/svm.peaks.$task.train.0.gc.seq $task/svm.peaks.$task.train.1.gc.seq $task/svm.peaks.$task.train.2.gc.seq $task/svm.peaks.$task.train.3.gc.seq $task/svm.peaks.$task.train.4.gc.seq $task/svm.peaks.$task.train.5.gc.seq $task/svm.peaks.$task.train.6.gc.seq $task/svm.peaks.$task.train.7.gc.seq $task/svm.peaks.$task.train.8.gc.seq $task/svm.peaks.$task.train.9.gc.seq
    
done
