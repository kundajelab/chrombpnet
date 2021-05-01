CUDA_VISIBLE_DEVICES=1 python score_model_overlap_merged.py --model_path /srv/scratch/annashch/chrombpnet/gm12878_dnase/with_bias_unplugged/gm12878.dnase.with.bias.unplugged.0.hdf5 \
		    --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.bpnet.unstranded.bw \
		    --out_prefix GM12878.DNASE.Fold0 \
		           --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz
