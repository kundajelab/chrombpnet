data_type=ATAC_PE

#uncorrected=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/data/GM12878/sorted_merged_uncorrected.bw 
#bias=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/data/GM12878/sorted_merged_expected.bw 
#corrected=/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/data/GM12878/sorted_merged_corrected.bw 
#orig=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/GM12878_unstranded.bw
uncorrected=HEPG2_ATAC_TOBIAS/sorted_merged_uncorrected.bw 
bias=HEPG2_ATAC_TOBIAS/sorted_merged_expected.bw 
corrected=HEPG2_ATAC_TOBIAS/sorted_merged_corrected.bw 
orig=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/HEPG2/data/HEPG2_unstranded.bw
outdir=HEPG2_max_$data_type

#motif=ETS1_HUMAN.H11MO.0.A_ETS1_HUMAN.H11MO.0.A
#motif=NRF1_MOUSE.H11MO.0.A_NRF1_MOUSE.H11MO.0.A
motif=HNF4G_MA0484.1_HNF4G_MA0484.1
#motif=FOXM1_HUMAN.H11MO.0.A_FOXM1_HUMAN.H11MO.0.A_1
#motif=ATF3_MOUSE.H11MO.0.A_ATF3_MOUSE.H11MO.0.A


#bound=$outdir/$motif/beds/$motif_HEPG2_bound.bed
#unbound=$outdir/$motif/beds/$motif_HEPG2_unbound.bed
bound=$outdir/$motif/beds/$motif"_HEPG2_bound.bed"
unbound=$outdir/$motif/beds/$motif"_HEPG2_unbound.bed"
motifname=HNF4G_try2


#TOBIAS FootprintScores --signal $corrected --regions $peaks --output $outdir/plots/tobias_fpd_scores.bw --cores 8


TOBIAS PlotAggregate --TFBS $outdir/$motif/beds/$motif"_all.bed" $outdir/$motif/beds/$motif"_HEPG2_bound.bed" $outdir/$motif/beds/$motif"_HEPG2_unbound.bed" --signals $uncorrected $bias $corrected $orig --output $outdir/plots/tobias_$motifname"_footprint.png" --share_y sites --plot_boundaries


mkdir $outdir/plots/

TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/tobias_ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6
TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/tobias_ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1



