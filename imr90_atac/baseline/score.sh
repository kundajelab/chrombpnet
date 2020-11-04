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

#for fold in `seq 0 4`
#do
#    kerasAC_score_bpnet_legacy \
#	--labels $outdir/$model_name.$fold.labels \
#	--predictions $outdir/$model_name.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf $outdir/$model_name.$fold.scores.smooth.labels \
#	--title "IMR90 ATAC, counts loss x 84, seed 1234" \
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
	--title "IMR90 ATAC, counts loss x 84, seed 1234" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/277549db-c2d8-49d3-ace0-81ad5d4088fb/call-align/IMR90.merged.pr1.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/277549db-c2d8-49d3-ace0-81ad5d4088fb/call-align/IMR90.merged.pr2.unstranded.bw
done
