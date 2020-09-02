for fold in `seq 0 9`
do
    kerasAC_score_bpnet \
	--labels predictions.H1ESC.withdups.20counts.1234seed.$fold.labels \
	--predictions predictions.H1ESC.withdups.20counts.1234seed.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.H1ESC.20counts.1234seed.$fold \
	--title "H1ESC, fold $fold, counts loss x20, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--
done
