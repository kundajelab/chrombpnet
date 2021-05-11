python predict_and_interpret_pwm_seqs.py --model_path /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
       --background_pickle GM12878_dnase.background.pickle \
       --pwm_consensus_pickle JASPER.top3.pickle \
       --motif_list_to_score NFKB1_MA0105.4 NFKB2_MA0778.1 Gata1_MA0035.3 GATA1+TAL1_MA0140.2 Hnf4a_MA0114.3 \
       --out_prefix gm12878.dnase.hit.list

       
