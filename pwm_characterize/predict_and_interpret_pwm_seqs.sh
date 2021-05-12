#python predict_and_interpret_pwm_seqs.py --model_path /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
#       --background_pickle GM12878_dnase.background.pickle \
#       --pwm_consensus_pickle JASPER.top3.pickle \
#       --motif_list_to_score NFKB1_MA0105.4 NFKB2_MA0778.1 Gata1_MA0035.3 GATA1+TAL1_MA0140.2 Hnf4a_MA0114.3 \
#       --out_prefix gm12878.dnase.hit.list


#CUDA_VISIBLE_DEVICES=0 python predict_and_interpret_pwm_seqs.py --model_path /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
#       --background_pickle K562_dnase.background.pickle \
#       --pwm_consensus_pickle JASPER.top3.pickle \
#       --motif_list_to_score NFKB1_MA0105.4 NFKB2_MA0778.1 Gata1_MA0035.3 GATA1+TAL1_MA0140.2 Hnf4a_MA0114.3 \
#       --out_prefix k562.dnase.hit.list

CUDA_VISIBLE_DEVICES=3 python predict_and_interpret_pwm_seqs.py --model_path /srv/scratch/annashch/chrombpnet/hepg2_dnase/with_bias_unplugged/hepg2.dnase.with.bias.unplugged.0.hdf5 \
       --background_pickle HEPG2_dnase.background.pickle \
       --pwm_consensus_pickle JASPER.top3.pickle \
       --motif_list_to_score NFKB1_MA0105.4 NFKB2_MA0778.1 Gata1_MA0035.3 GATA1+TAL1_MA0140.2 Hnf4a_MA0114.3 \
       --out_prefix hepg2.dnase.hit.list
