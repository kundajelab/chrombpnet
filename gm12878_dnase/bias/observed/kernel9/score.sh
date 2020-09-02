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
	--title "GM12878, fold $fold, counts loss x64, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/atac/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.pr2.bam.bpnet.unstranded.bw
done

