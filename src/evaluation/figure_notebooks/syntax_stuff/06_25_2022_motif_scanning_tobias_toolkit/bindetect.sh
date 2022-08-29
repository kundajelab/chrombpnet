data_type=ATAC_PE
#data_type=DNASE_PE
model_name=HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
#model_name=HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0

#mode=counts
mode=profile

#motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/06_22_2022_motif_scanning/ppms/$mode.ppms.txt
#motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/ppms/$mode.ppms.txt
motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/ppms/counts.ppms.txt

model_name=nautilus_runs_jun16/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
signal=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/merged.HEPG2.$mode.bw

uncorrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2_w_bias.bw
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2bias.bw
corrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2_wo_bias.bw
outdir=output/HEPG2_mean_$mode"_"$data_type

# same for both atac and dnase

peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/merged.HEPG2.interpreted_regions.bed
awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' $peaks > temp.bed
peaks=temp.bed
peaks_header=header.txt
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa

#TOBIAS BINDetect --motifs $motifs --signals $signal --genome $genome --peaks $peaks --peak_header $peaks_header --outdir $outdir --cond_names HEPG2 --cores 8

mkdir  $outdir/plots/

#motif=ETS1_HUMAN.H11MO.0.A_ETS1_HUMAN.H11MO.0.A
#motif=NRF1_MOUSE.H11MO.0.A_NRF1_MOUSE.H11MO.0.A
motif=NRF1_HUMAN.H11MO.0.A_NRF1_HUMAN.H11MO.0.A

#motif=HNF4G_MA0484.1_HNF4G_MA0484.1
#motif=FOXM1_HUMAN.H11MO.0.A_FOXM1_HUMAN.H11MO.0.A_1
#motif=ATF3_MOUSE.H11MO.0.A_ATF3_MOUSE.H11MO.0.A

#bound=$outdir/$motif/beds/$motif_HEPG2_bound.bed
#unbound=$outdir/$motif/beds/$motif_HEPG2_unbound.bed
bound=$outdir/$motif/beds/$motif"_HEPG2_bound.bed"
unbound=$outdir/$motif/beds/$motif"_HEPG2_unbound.bed"

motifname=$1
motif=$motifname"_"$motifname


#motifname=FOXM1
#motifname=ATF3
#motifname=HEPG2
#motifname=NRF1

#TOBIAS FootprintScores --signal $corrected --regions $peaks --output $outdir/plots/fpd_scores.bw --cores 8

#TOBIAS PlotAggregate --TFBS $outdir/$motif/beds/$motif"_all.bed" $outdir/$motif/beds/$motif"_HEPG2_bound.bed" $outdir/$motif/beds/$motif"_HEPG2_unbound.bed" --signals $uncorrected $bias $corrected $orig --output $outdir/plots/$motifname"_footprint.png" --share_y sites --plot_boundaries


#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6
#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1

TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1
TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6










