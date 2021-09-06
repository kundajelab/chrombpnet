outdir=$1
model_name=$2
fold=$3
cell_line=$4
seed=$5
kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.inpeak.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.inpeak.scores \
    --title "$cell_line only bias, fold $fold, seed $seed" \
    --label_min_to_score 2.3 \
    --label_max_to_score 11.5 \
    --num_tasks 1

kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.inpeak.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.inpeak.scores.smooth.labels \
    --title "$cell_line only bias, fold $fold, seed $seed" \
    --label_min_to_score 2.3 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps


kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.inpeak.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.inpeak.scores.smooth.both \
    --title "$cell_line only bias, fold $fold, seed $seed" \
    --label_min_to_score 2.3 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps \
    --smooth_predicted_profile

