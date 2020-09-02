kerasAC_score_bpnet \
    --labels predictions.k562.withdups.1234seed.250counts.500filters.0.labels \
    --predictions predictions.k562.withdups.1234seed.250counts.500filters.0.predictions \
    --losses profile counts \
    --loss_suffixes 0 1 \
    --outf score.k562.250counts.1234seed.0.500filters \
    --title "K562, fold 0, counts loss x250, seed 1234, 500 filters"
