#data_type=ATAC_PE
data_type=DNASE_PE
#model_name=K562_02.17.2022_bias_128_4_1234_0.5_fold_0
mode=counts

motifs=/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/06_22_2022_motif_scanning/ppms/counts.ppms.txt

model_name=nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0

signal=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/K562/$model_name/interpret/merged.K562.$mode.bw
uncorrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/K562/$model_name/interpret/K562_w_bias.bw
bias=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/K562/$model_name/interpret/K562bias.bw
corrected=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/K562/$model_name/interpret/K562_wo_bias.bw

outdir=output/K562_mean_$mode'_'$data_type

# same for both atac and dnase

peaks=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/$data_type/K562/$model_name/interpret/merged.K562.interpreted_regions.bed
awk ' {print $1 "\t" $2 + $10 - 500 "\t" $2 + $10 + 500}' $peaks > k562_temp.bed
peaks=k562_temp.bed
peaks_header=header.txt
genome=/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa

TOBIAS BINDetect --motifs $motifs --signals $signal --genome $genome --peaks $peaks --peak_header $peaks_header --outdir $outdir --cond_names K562 --cores 8


motifname=$1
motif=$motifname"_"$motifname


mkdir  $outdir/plots/


#TOBIAS PlotAggregate --TFBS $outdir/$motif/beds/$motif"_all.bed" $outdir/$motif/beds/$motif"_K562_bound.bed" $outdir/$motif/beds/$motif"_K562_unbound.bed" --signals $uncorrected $bias $corrected --output $outdir/plots/NRF1_footprint.png --share_y sites --plot_boundaries
bound=$outdir/$motif/beds/$motif"_K562_bound.bed"
unbound=$outdir/$motif/beds/$motif"_K562_unbound.bed"

#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/ranked_by_affinity_NRF1_heatmap.png --signal_labels NRF1 --share_colorbar --sort_by -6
#TOBIAS PlotHeatmap --TFBS $bound $unbound --signals $corrected --output $outdir/plots/ranked_by_contrib_NRF1_heatmap.png --signal_labels NRF1 --share_colorbar --sort_by -1

#TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/ranked_by_contrib_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -1
#TOBIAS PlotHeatmap --TFBS  $outdir/$motif/beds/$motif"_all.bed" --signals $corrected --output $outdir/plots/ranked_by_affinity_$motifname"_heatmap.png" --signal_labels HNF4G --share_colorbar --sort_by -6





