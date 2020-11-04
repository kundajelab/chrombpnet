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
#	--title "H1ESC ATAC, counts loss x 56, seed 1234" \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5
#done

#for fold in `seq 0 4`
#do
#    kerasAC_score_bpnet_legacy \
#	--labels $outdir/$model_name.$fold.labels \
#	--predictions $outdir/$model_name.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf $outdir/$model_name.$fold.scores.smooth.labels \
#	--title "H1ESC ATAC, counts loss x 56, seed 1234" \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5 \
#	--smooth_observed_profile \
#	--smooth_preps
#done
for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels $outdir/$model_name.$fold.labels \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf $outdir/$model_name.$fold.scores.smooth.both \
	--title "H1ESC ATAC, counts loss x 56, seed 1234" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/58fb3f13-be45-45de-8a39-d0bfbeaf86c5/call-align/H1ESC.merged.pr1.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/58fb3f13-be45-45de-8a39-d0bfbeaf86c5/call-align/H1ESC.merged.pr2.unstranded.bw
done
