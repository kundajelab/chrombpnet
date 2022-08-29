data_type=ATAC_PE
model_name=HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
#model_name=HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0

#motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/06_22_2022_motif_scanning/ppms/counts.ppms.txt
motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/06_22_2022_motif_scanning/ppms/counts.ppms.txt

model_name=nautilus_runs_jun16/HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0
signal=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/merged.HEPG2.counts.bw
uncorrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2_w_bias.bw
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2bias.bw
corrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/HEPG2_wo_bias.bw
outdir=HEPG2_max_$data_type

# same for both atac and dnase

peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/HEPG2/$model_name/interpret/merged.HEPG2.interpreted_regions.bed
awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' $peaks > temp.bed
peaks=temp.bed
peaks_header=header.txt
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa

#TOBIAS BINDetect --motifs $motifs --signals $signal --genome $genome --peaks $peaks --peak_header $peaks_header --outdir $outdir --cond_names HEPG2 --cores 8


mkdir  $outdir/plots/

TOBIAS PlotAggregate --TFBS HEPG2_max_DNASE_PE/HNF4G_MA0484.1_HNF4G_MA0484.1/beds/HNF4G_MA0484.1_HNF4G_MA0484.1_all.bed HEPG2_max_DNASE_PE/HNF4G_MA0484.1_HNF4G_MA0484.1/beds/HNF4G_MA0484.1_HNF4G_MA0484.1_HEPG2_bound.bed HEPG2_max_DNASE_PE/HNF4G_MA0484.1_HNF4G_MA0484.1/beds/HNF4G_MA0484.1_HNF4G_MA0484.1_HEPG2_unbound.bed --signals $uncorrected $bias $corrected --output $outdir/plots/dnase_labels_HNF4G_footprint.png --share_y sites --plot_boundaries
bound=HEPG2_max_DNASE_PE/HNF4G_MA0484.1_HNF4G_MA0484.1/beds/HNF4G_MA0484.1_HNF4G_MA0484.1_HEPG2_bound.bed
unbound=HEPG2_max_DNASE_PE/HNF4G_MA0484.1_HNF4G_MA0484.1/beds/HNF4G_MA0484.1_HNF4G_MA0484.1_HEPG2_unbound.bed

TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/dnase_labels_ranked_by_affinity_HNF4G_heatmap.png --signal_labels HNF4G --share_colorbar --sort_by -6
TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/dnase_labels_ranked_by_contrib_HNF4G_heatmap.png --signal_labels HNF4G --share_colorbar --sort_by -1
#TOBIAS FootprintScores --signal $corrected --regions $peaks --output $outdir/plots/fpd_scores.bw --cores 8





