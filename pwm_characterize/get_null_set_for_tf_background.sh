#GM12878 DNASE, use overlap peaks for background
peak_file="/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz"
python get_null_set_for_tf_background.py \
       --peak_file $peak_file \
       --n_to_sample 100 \
       --out_pickle GM12878_dnase.background.pickle \
       --bigwig_count_file /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.bpnet.unstranded.bw \
       --model_path /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
 
