outdir=$1
model_name=$2
fold=$3
cell_line=$4
seed=$5
min_logcount=$6
max_logcount=$7

kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores \
    --title "$cell_line bias finetuning, fold $fold, seed $seed" \
    --label_min_to_score $min_logcount \
    --label_max_to_score $max_logcount \
    --num_tasks 1

kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.labels \
    --title "$cell_line bias finetuning, fold $fold, seed $seed" \
    --label_min_to_score $min_logcount \
    --label_max_to_score $max_logcount \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps


kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.both \
    --title "$cell_line bias finetuning, fold $fold, seed $seed" \
    --label_min_to_score $min_logcount \
    --label_max_to_score $max_logcount \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps \
    --smooth_predicted_profile

