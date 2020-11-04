outdir=$1
model_name=$2

#for fold in `seq 0 4`
#do
#    kerasAC_score_bpnet_legacy \
#	--labels $outdir/$model_name.$fold.labels \
#	--predictions $outdir/$model_name.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf $outdir/$model_name.$fold.scores.smoothed.labels \
#	--title "GM12878 ATAC, counts loss x 55, seed 1234" \
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
	--outf $outdir/$model_name.$fold.scores.smoothed.both \
	--title "GM12878 ATAC, counts loss x 55, seed 1234" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-align/GM12878.merged.pr1.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-align/GM12878.merged.pr2.unstranded.bw
done
