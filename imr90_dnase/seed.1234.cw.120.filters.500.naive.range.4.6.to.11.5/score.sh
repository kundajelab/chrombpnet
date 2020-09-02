#for fold in 0 `seq 0 4`
#do 
#    kerasAC_score_bpnet_legacy \
#	--labels seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.labels \
#	--predictions seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.predictions \
#	--losses profile counts \
#	--loss_suffixes 0 1 \
#	--outf score.imr90.101.1234seed.$fold \
#	--title "IMR90, fold $fold, counts loss x120, seed 1234" \
#	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr1.sorted.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr2.bam.bpnet.unstranded.bw \
#	--label_min_to_score 4.6 \
#	--label_max_to_score 11.5
#done

for fold in 0 `seq 0 4`
do 
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.labels \
	--predictions seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.imr90.101.1234seed.$fold.smooth.labels \
	--title "IMR90, fold $fold, counts loss x120, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr1.sorted.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_preps
done

for fold in 0 `seq 0 4`
do 
    kerasAC_score_bpnet_legacy \
	--labels seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.labels \
	--predictions seed.1234.cw.120.filters.500.naive.range.4.6.to.11.5.$fold.predictions \
	--losses profile counts \
	--loss_suffixes 0 1 \
	--outf score.imr90.101.1234seed.$fold.smooth.both \
	--title "IMR90, fold $fold, counts loss x120, seed 1234" \
	--pseudoreps /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr1.sorted.bam.bpnet.unstranded.bw /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/dnase/f38bfd43-b57f-4c55-be06-b02d3f16512a/call-bowtie2/shard-0/execution/ENCSR477RTP.merged.bam.pr2.bam.bpnet.unstranded.bw \
	--label_min_to_score 4.6 \
	--label_max_to_score 11.5 \
	--smooth_observed_profile \
	--smooth_predicted_profile \
	--smooth_preps
done

