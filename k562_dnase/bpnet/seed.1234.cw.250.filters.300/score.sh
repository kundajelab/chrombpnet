for fold in `seq 0 9`
do
    kerasAC_score_bpnet \
	--labels predictions.k562.withdups.1234seed.250counts.$fold.labels \
	--predictions predictions.k562.withdups.1234seed.250counts.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.k562.250counts.1234seed.$fold \
	--title "K562, fold $fold, counts loss x250, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr2.bam.bpnet.unstranded.bw    
done
