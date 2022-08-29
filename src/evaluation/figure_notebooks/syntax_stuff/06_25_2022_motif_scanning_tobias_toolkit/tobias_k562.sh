input_bam=/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/K562/sorted_merged.bam
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa
peaks=k562_temp.bed
blacklist=/mnt/lab_data2/anusri/chrombpnet/reference/blacklist_slop1057.bed
out_dir=K562_ATAC_TOBIAS

#TOBIAS ATACorrect --bam $input_bam --genome $genome --peaks $peaks --blacklist $blacklist --outdir $out_dir --cores 8

corrected=K562_ATAC_TOBIAS/sorted_merged_corrected.bw 
uncorrected=K562_ATAC_TOBIAS/sorted_merged_uncorrected.bw 
expected=K562_ATAC_TOBIAS/sorted_merged_expected.bw 

peaks=k562_temp.bed

#TOBIAS FootprintScores --signal $corrected --regions $peaks --output $out_dir/tobias_fpd_scores.bw --cores 8


motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/06_22_2022_motif_scanning/ppms/counts.ppms.txt
signal=$out_dir/tobias_fpd_scores.bw
peaks_header=header.txt

TOBIAS BINDetect --motifs $motifs --signals $signal --genome $genome --peaks $peaks --peak_header $peaks_header --outdir $out_dir --cond_names K562 --cores 8

#motif=ETS1_HUMAN.H11MO.0.A_ETS1_HUMAN.H11MO.0.A
#motif=NRF1_MOUSE.H11MO.0.A_NRF1_MOUSE.H11MO.0.A
#motif=NRF1_HUMAN.H11MO.0.A_NRF1_HUMAN.H11MO.0.A
#motif=HNF4G_MA0484.1_HNF4G_MA0484.1
#motif=FOXM1_HUMAN.H11MO.0.A_FOXM1_HUMAN.H11MO.0.A_1
#motif=ATF3_MOUSE.H11MO.0.A_ATF3_MOUSE.H11MO.0.A

#bound=$outdir/$motif/beds/$motif_K562_bound.bed
#unbound=$outdir/$motif/beds/$motif_K562_unbound.bed
outdir=$out_dir

bound=$outdir/$motif/beds/$motif"_K562_bound.bed"
unbound=$outdir/$motif/beds/$motif"_K562_unbound.bed"
#motifname=HNF4G
#motifname=NRF1
motifname=$1
motif=$motifname"_"$motifname




#TOBIAS FootprintScores --signal $corrected --regions $peaks --output $outdir/plots/tobias_fpd_scores.bw --cores 8

mkdir $outdir/plots/

#TOBIAS PlotAggregate --TFBS $outdir/$motif/beds/$motif"_all.bed" $outdir/$motif/beds/$motif"_K562_bound.bed" $outdir/$motif/beds/$motif"_K562_unbound.bed" --signals $uncorrected $bias $corrected $orig --output $outdir/plots/tobias_$motifname"_footprint.png" --share_y sites --plot_boundaries

#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/tobias_ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6
#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/tobias_ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1

#TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/tobias_ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1
#TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/tobias_ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6



