#!/bin/bash
./train.sh 0 1 6mer.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/seed2345 params.6mer.txt
#./train.counts.sh 0 1 counts.6mer.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.6mer.txt
#./train.profile.sh 0 1 profile.6mer.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.6mer.txt
#!/bin/bash
./train.sh 0 1 bpnet.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/seed2345 params.bpnet.txt
#./train.counts.sh 0 1 counts.bpnet.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.bpnet.txt
#./train.profile.sh 0 1 profile.bpnet.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.bpnet.txt
#!/bin/bash
./train.sh 0 1 tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias/seed2345 params.tobias.txt
#./train.counts.sh 0 1 counts.tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.tobias.txt
#./train.profile.sh 0 1 profile.tobias.bias.preds.hepg2 2345 /srv/scratch/annashch/deeplearning/profile/hepg2_atac/pred_from_bias params.tobias.txt
