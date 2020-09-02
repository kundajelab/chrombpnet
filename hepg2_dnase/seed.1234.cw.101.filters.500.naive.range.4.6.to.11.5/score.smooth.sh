for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cs.101.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.labels \
	--predictions seed.1234.cs.101.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.smooth.labels.hepg2.101.1234seed.$fold \
	--title "HEPG2, fold $fold, counts loss x101, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_preps
done

for fold in `seq 0 4`
do
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cs.101.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.labels \
	--predictions seed.1234.cs.101.filters.500.naive.range.4.6.to.11.5.bias.observed.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.smooth.both.hepg2.101.1234seed.$fold \
	--title "HEPG2, fold $fold, counts loss x101, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.pr1.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/38f0a76b-e6c6-444e-84e5-b5a98a554694/call-bowtie2/shard-0/execution/ENCSR149XIL.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_preps \
	--smooth_predicted_profile
done

