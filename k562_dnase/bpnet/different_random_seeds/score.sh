fold=0
for seed in 1234 #2345 3456 4567 5678
do
    kerasAC_score_bpnet \
	--labels predictions.k562.withdups.$seed\seed.300counts.$fold.labels \
	--predictions predictions.k562.withdups.$seed\seed.300counts.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.k562.300counts.$seed\seed.$fold \
	--title "K562, fold $fold, counts loss x300, seed $seed"
done
