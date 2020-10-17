kerasAC_score_bpnet \
    --predictions preds.hdf5.predictions \
    --outf scores \
    --title "K562 H3K27ac, counts loss x 25, seed 1234, lr=0.005" \
    --label_min_to_score 3 \
    --label_max_to_score 11.5 \
    --num_tasks 2 \
    --losses profile counts \
    --smooth_observed_profile \
    --smooth_predicted_profile


