#GM12878
CUDA_VISIBLE_DEVICES=1 python score_model_overlap_merged.py --model_path /srv/scratch/annashch/chrombpnet/gm12878_atac/with_bias_unplugged/gm12878.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/5846e593-a935-4bd9-9294-422a05f9f9b8/call-align/GM12878.merged.bam.bpnet.unstranded.bw \
       --out_prefix GM12878.ATAC.Overlap.Fold0 \
       --bed_file_to_score /srv/scratch/annashch/chrombpnet/score/overlap_peaks/GM12878.overlap.ATAC.DNASE.merged.formatted.for.scoring.bed \
       --precentered_intervals \
       --presorted_intervals

#HEPG2
CUDA_VISIBLE_DEVICES=1 python score_model_overlap_merged.py --model_path /srv/scratch/annashch/chrombpnet/hepg2_atac/with_bias_unplugged/hepg2.atac.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/atac/caper_out/25b3429e-5864-4e8d-a475-a92df8938887/call-align/HEPG2.merged.bam.bpnet.unstranded.bw \
       --out_prefix HEPG2.ATAC.Overlap.Fold0 \
       --bed_file_to_score /srv/scratch/annashch/chrombpnet/score/overlap_peaks/HEPG2.overlap.ATAC.DNASE.merged.formatted.for.scoring.bed \
       --precentered_intervals \
       --presorted_intervals
