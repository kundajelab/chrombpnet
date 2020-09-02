outdir=$1
model_name=$2
fold=$3
kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores \
    --title "HEPG2 DNASE with BPNET bias, fold $fold,  seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1
kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.labels \
    --title "HEPG2 DNASE with BPNET bias, fold $fold,  seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps

kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.both \
    --title "HEPG2 DNASE with BPNET bias, fold $fold,  seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_predicted_profile \
    --smooth_preps 

