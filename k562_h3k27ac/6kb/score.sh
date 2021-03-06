outdir=$1
model_name=$2

for fold in 0 #`seq 0 4`
do
    kerasAC_score_bpnet \
	--predictions $outdir/$model_name.$fold.predictions \
	--outf $outdir/$model_name.$fold.scores \
	--title "K562 H3K27ac, counts loss x 25, seed 1234" \
	--label_min_to_score 3 \
	--label_max_to_score 11.5 \
	--num_tasks 2 \
	--losses profile counts
done
