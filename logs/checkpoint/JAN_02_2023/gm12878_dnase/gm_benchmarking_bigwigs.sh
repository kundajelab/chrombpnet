regions1=GM12878.interpreted_regions.bed

out=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/gm12878_dnase/

model=/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR000EMT/chrombpnet_model_feb15_fold_0/chrombpnet_wo_bias.h5
#celline=gmcorrected

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

#mkdir $output_dir
#bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 


#model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/HEPG2/copy_HEPG2_06.08.2022_bias_128_4_1234_0.8_fold_0/bias_model/bias.h5
#celline=hepbias

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
#bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 



model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5
celline=corrected

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 


model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/bias_model/bias.h5
celline=bias

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 


model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/uncorrected_model_05.10.2022/uncorrected_model/hint_atac.h5 
celline=uncorrected

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 


model=/mnt/lab_data2/anusri/chrombpnet/results/hint_atac/DNASE_SE/GM12878/DNASE_SE_11.28.2022_hint_atac/hint_atac_model/hint_atac.h5
celline=hint

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 


model=/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/DNASE_SE_03.06.2022_simplebias/chrombpnet_model/chrombpnet_wo_bias.h5
celline=simplebias

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=2
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 




