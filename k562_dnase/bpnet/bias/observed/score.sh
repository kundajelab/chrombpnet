outdir=$1
model_name=$2

for fold in `seq 0 4`
do
    kerasAC_score_bpnet \
	--labels $outdir/$model_name.$fold.labels \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf $outdir/$model_name.$fold.scores \
	--title "K562, fold $fold, counts loss x250, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr2.bam.bpnet.unstranded.bw    
done
