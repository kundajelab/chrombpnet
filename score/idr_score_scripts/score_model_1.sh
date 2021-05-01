#GM12878
CUDA_VISIBLE_DEVICES=1 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/gm12878_atac/with_bias_unplugged/gm12878.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-align/GM12878.merged.bam.bpnet.unstranded.bw \
       --out_prefix GM12878.ATAC.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
#HEPG2
CUDA_VISIBLE_DEVICES=1 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/hepg2_atac/with_bias_unplugged/hepg2.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-align/HEPG2.merged.bam.bpnet.unstranded.bw \
       --out_prefix HEPG2.ATAC.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-reproducibility_idr/glob-1b1244d5baf1a7d98d4b7b76d79e43bf/idr.optimal_peak.narrowPeak.gz
