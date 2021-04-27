#K562
CUDA_VISIBLE_DEVICES=0 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/with_bias_unplugged/k562.dnase.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.bpnet.unstranded.bw \
       --out_prefix K562.DNASE.Fold0 \
       --bed_file_to_score /srv/scratch/annashch/chrombpnet/k562_dnase/bpnet/K562.dnase.idr.optimal_peak.narrowPeak.gz
#HEPG2
CUDA_VISIBLE_DEVICES=0 python score_model.py --model_path /srv/scratch/annashch/chrombpnet/hepg2_dnase/with_bias_unplugged/hepg2.dnase.with.bias.unplugged.0.hdf5 \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.bpnet.unstranded.bw \
       --out_prefix HEPG2.DNASE.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/ce805260-55f8-43c8-b2a1-a232b4a0e369/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz

#GM12878
CUDA_VISIBLE_DEVICES=0 python score_model.py --model_path /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.bpnet.unstranded.bw \
       --bigwig_labels /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-bowtie2/shard-0/execution/ENCSR000EMT.merged.bam.bpnet.unstranded.bw \
       --out_prefix GM12878.DNASE.Fold0 \
       --bed_file_to_score /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/13da5ebe-0941-4855-8599-40bbcc5c58b4/call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz


