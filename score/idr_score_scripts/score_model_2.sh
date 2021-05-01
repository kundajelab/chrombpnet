#H1HESC
CUDA_VISIBLE_DEVICES=2 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/h1esc_atac/with_bias_unplugged/h1esc.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/58fb3f13-be45-45de-8a39-d0bfbeaf86c5/call-align/H1ESC.merged.bam.bpnet.unstranded.bw \
       --out_prefix H1HESC.ATAC.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/58fb3f13-be45-45de-8a39-d0bfbeaf86c5/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz

#IMR90
CUDA_VISIBLE_DEVICES=2 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/imr90_atac/with_bias_unplugged/imr90.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/277549db-c2d8-49d3-ace0-81ad5d4088fb/call-align/IMR90.merged.bam.bpnet.unstranded.bw \
       --out_prefix IMR90.ATAC.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/277549db-c2d8-49d3-ace0-81ad5d4088fb/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
