outdir=$1
model_name=$2
fold=$3
kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores \
    --title "DNASE with bias, fold $fold, counts loss x18, seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr1.bam.bpnet.unstranded.bw,/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr2.bam.bpnet.unstranded.bw

kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.labels \
    --title "DNASE with bias, fold $fold, counts loss x18, seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_preps \
    --pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr1.bam.bpnet.unstranded.bw,/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr2.bam.bpnet.unstranded.bw


kerasAC_score_bpnet \
    --predictions $outdir/$model_name.$fold.predictions \
    --losses profile counts \
    --outf $outdir/$model_name.$fold.scores.smooth.both \
    --title "DNASE with bias, fold $fold, counts loss x18, seed 1234" \
    --label_min_to_score 4.6 \
    --label_max_to_score 11.5 \
    --num_tasks 1 \
    --smooth_observed_profile \
    --smooth_predicted_profile \
    --smooth_preps \
    --pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr1.bam.bpnet.unstranded.bw,/oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/09ce5f39-5360-411b-88dd-b86f4a1286a7/call-bowtie2/shard-0/execution/ENCSR000EOT.merged.bam.pr2.bam.bpnet.unstranded.bw

