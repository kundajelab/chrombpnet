CUDA_VISIBLE_DEVICES=0 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.bpnet.unstranded.bw \
       --out_prefix K562.DNASE.Fold0 \
       --bed_file_to_score /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/K562.dnase.idr.optimal_peak.narrowPeak.gz


