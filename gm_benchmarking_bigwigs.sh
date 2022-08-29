chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/bias_model_scaled.h5
celline=GM12878
gpu=1
regions1=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/interpret/merged.GM12878.interpreted_regions.bed

output_dir=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/interpret/

mkdir $output_dir

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir

chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
bias=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
celline=GM12878
gpu=1
regions1=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/interpret/merged.GM12878.interpreted_regions.bed

output_dir=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/interpret/

mkdir $output_dir

bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir




regions=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/interpret/GM12878.interpreted_regions.merged.bed
bigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/merged/GM12878/ATAC/GM12878_unstranded.bw
obigwig=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/GM12878_unstranded_scaled.bw

#python normalize_bigwigs_from_peaks.py --input_bed $regions -bigwig $bigwig --output_path $obigwig


chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
bias=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5
celline=uncorrected
gpu=1
#regions=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/interpret/GM12878.interpreted_regions.bed
regions=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/interpret/GM12878.interpreted_regions.merged.bed

output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/uncorrected/
interpret_dir=results/chrombpnet/ATAC_PE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/interpret/GM12878

mkdir $output_dir
#bash make_bigwig_new_2.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir $interpret_dir


chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/chrombpnet.h5
celline=simplebias
gpu=1
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/simplebias/
interpret_dir=results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_03.06.2022_simplebias/chrombpnet_model/interpret/GM12878

mkdir $output_dir
#bash make_bigwig_new_2.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir $interpret_dir

chrombpnet_nb=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
chrombpnet=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet.h5
bias=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/bias_model_scaled.h5
celline=chrombpnet
gpu=1
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/chrombpnet/
interpret_dir=results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/interpret/GM12878
mkdir $output_dir
#bash make_bigwig_new_2.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir $interpret_dir

chrombpnet_nb=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/hint_atac.h5
chrombpnet=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/hint_atac.h5
bias=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/hint_atac.h5
celline=hint
gpu=1
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/hint/
interpret_dir=results/hint_atac/ATAC_PE/GM12878/ATAC_PE_03.06.2022_hint_atac/hint_atac_model/interpret/GM12878

mkdir $output_dir
#bash make_bigwig_new_2.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir $interpret_dir 

chrombpnet_nb=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias_wo_bias.h5
chrombpnet=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias_wo_bias.h5
bias=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/tobias_wo_bias.h5
celline=tobias
gpu=1
output_dir=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/benchmarking/tobias/
interpret_dir=results/tobias/ATAC_PE/GM12878/ATAC_PE_04.06.2022_tobias_with_bias_bigwig/tobias_model/interpret/GM12878

mkdir $output_dir
#bash make_bigwig_new_2.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions $output_dir $interpret_dir
