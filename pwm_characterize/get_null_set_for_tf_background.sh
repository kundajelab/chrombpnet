#GM12878 DNASE, use overlap peaks for background
#peak_file="/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz"
#python get_null_set_for_tf_background.py \
#       --peak_file $peak_file \
#       --n_to_sample 100 \
#       --out_pickle GM12878_dnase.background.pickle \
#       --bigwig_count_file /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.bpnet.unstranded.bw \
#       --model_path /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
 
#K562 DNASE
#peak_file=/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz
#CUDA_VISIBLE_DEVICES=0 python get_null_set_for_tf_background.py \
#       --peak_file $peak_file \
#       --n_to_sample 100 \
#       --out_pickle K562_dnase.background.pickle \
#       --model_path /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
#       --bigwig_count_file /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.bpnet.unstranded.bw &

#HEPG2 DNASE
peak_file=/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/ce805260-55f8-43c8-b2a1-a232b4a0e369/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz
CUDA_VISIBLE_DEVICES=3 python get_null_set_for_tf_background.py \
       --peak_file $peak_file \
       --n_to_sample 100 \
       --out_pickle HEPG2_dnase.background.pickle \
       --model_path /srv/scratch/annashch/chrombpnet/hepg2_dnase/with_bias_unplugged/hepg2.dnase.with.bias.unplugged.0.hdf5 \
       --bigwig_count_file /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.bpnet.unstranded.bw 
