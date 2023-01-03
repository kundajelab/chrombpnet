output_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/footprints_motifs_uncorrected/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=chrombpnet.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
#motifs=AP1,CTCF
model_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/
mkdir $output_dir

gpu=1

CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir

output_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/bias_model/footprints_motifs/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=bias.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
model_dir=results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/bias_model/

mkdir $output_dir
CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir



output_dir=results/chrombpnet/DNASE_SE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/footprints_motifs/
reference_fasta=reference/hg38.genome.fa
fold=splits/fold_0.json 
model=hint_atac.h5
motifs=tn5_1,tn5_2,tn5_3,tn5_4,tn5_5,GATA+TAL,AP1,CTCF,ETS,RUNX,NRF1,NFKB,SPI1
model_dir=results/chrombpnet/DNASE_SE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/

mkdir $output_dir
CUDA_VISIBLE_DEVICES=$gpu bash marginal_footprinting_all.sh $output_dir $reference_fasta $fold $model $motifs $model_dir



