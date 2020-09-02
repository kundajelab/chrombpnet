outdir=$1
model_name=$2

#for fold in 0 
#do
#    kerasAC_score_bpnet \
#	--predictions $outdir/$model_name.$fold.predictions \
#	--losses profile counts \
#	--outf $outdir/$model_name.$fold.scores \
#	--title "HEPG2 from $2" \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5 \
#	--num_tasks 1
#done

for fold in 0 
do
    kerasAC_score_bpnet \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--outf $outdir/$model_name.$fold.scores.smooth.labels \
	--title "HEPG2 from $2" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--num_tasks 1 \
	--smooth_observed_profile \
	--smooth_preps
done

for fold in 0 
do
    kerasAC_score_bpnet \
	--predictions $outdir/$model_name.$fold.predictions \
	--losses profile counts \
	--outf $outdir/$model_name.$fold.scores.smooth.both \
	--title "HEPG2 from $2" \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--num_tasks 1 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps
done
