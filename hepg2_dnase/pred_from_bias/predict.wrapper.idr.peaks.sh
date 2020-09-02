#!/bin/bash
#SEED 1234 
#./predict.idr.peaks.sh 0 0 seed1234/bpnet.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds &
#./predict.idr.peaks.sh 0 1 seed1234/tobias.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds &
#./predict.idr.peaks.sh 0 2 seed1234/6mer.bias.preds.hepg2 1234 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds &

#SEED 2345 
./predict.idr.peaks.sh 0 0 seed2345/bpnet.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 seed2345/tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 seed2345/6mer.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds


#SEED 3456
./predict.idr.peaks.sh 0 0 seed3456/bpnet.bias.preds.hepg2 3456 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 seed3456/tobias.bias.preds.hepg2 3456 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds
./predict.idr.peaks.sh 0 0 seed3456/6mer.bias.preds.hepg2 3456 /srv/scratch/annashch/deeplearning/profile/hepg2_dnase/pred_from_bias/idr_preds
