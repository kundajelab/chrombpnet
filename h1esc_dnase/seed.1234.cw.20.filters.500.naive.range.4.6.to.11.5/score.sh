#for fold in `seq 0 4`
#do
#    kerasAC_score_bpnet_legacy \
#	--labels seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.labels \
#	--predictions seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf score.h1esc.20counts.1234seed.$fold \
#	--title "H1ESC, fold $fold, counts loss x20, seed 1234" \
#	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr2.bam.bpnet.unstranded.bw \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5
#done
for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.labels \
	--predictions seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.h1esc.20counts.1234seed.$fold.smooth.labels \
	--title "H1ESC, fold $fold, counts loss x20, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_preps 
done
for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.labels \
	--predictions seed.1234.cs.20.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.h1esc.20counts.1234seed.$fold.smooth.both \
	--title "H1ESC, fold $fold, counts loss x20, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/c9ef8473-1374-41ef-9fab-8f07288e94e7/call-bowtie2/shard-0/execution/ENCSR000EMU.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps
done
