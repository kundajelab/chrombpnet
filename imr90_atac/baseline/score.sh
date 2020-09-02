outdir=$1
model_name=$2

#for fold in `seq 0 4`
#do
#    kerasAC_score_bpnet_legacy \
#	--labels $outdir/$model_name.$fold.labels \
#	--predictions $outdir/$model_name.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf $outdir/$model_name.$fold.scores \
#	--title "IMR90 ATAC, counts loss x 84, seed 1234" \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5
#done

for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels $outdir/$model_name.$fold.labels \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf $outdir/$model_name.$fold.scores.smooth.labels \
	--title "IMR90 ATAC, counts loss x 84, seed 1234" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_preps
done
for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels $outdir/$model_name.$fold.labels \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf $outdir/$model_name.$fold.scores.smooth.both \
	--title "IMR90 ATAC, counts loss x 84, seed 1234" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps
done
