#H1HESC
CUDA_VISIBLE_DEVICES=3 python score_model_overlap_merged.py --model_path /srv/scratch/annashch/chrombpnet/h1esc_dnase/with_bias_unplugged/h1esc.dnase.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.bpnet.unstranded.bw \
       --out_prefix H1HESC.DNASE.Overlap.Fold0 \
       --bed_file_to_score /srv/scratch/annashch/chrombpnet/score/overlap_peaks/H1HESC.overlap.ATAC.DNASE.merged.formatted.for.scoring.bed \
       --precentered_intervals \
       --presorted_intervals

#IMR90
CUDA_VISIBLE_DEVICES=3 python score_model_overlap_merged.py --model_path /srv/scratch/annashch/chrombpnet/imr90_dnase/with_bias_unplugged/imr90.dnase.with.bias.unplugged.0.hdf5 \
		    --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.bpnet.unstranded.bw \
		    --out_prefix IMR90.DNASE.Overlap.Fold0 \
		    --bed_file_to_score /srv/scratch/annashch/chrombpnet/score/overlap_peaks/IMR90.overlap.ATAC.DNASE.merged.formatted.for.scoring.bed \
		    --precentered_intervals \
		    --presorted_intervals
