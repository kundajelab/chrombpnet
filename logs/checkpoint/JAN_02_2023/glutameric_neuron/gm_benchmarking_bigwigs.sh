regions1=regions.bed

out=/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/glutameric_neuron/

model=/oak/stanford/groups/akundaje/salil512/chd/model_training/chrombpnet/fetal_brain/chrombpnet_models/c1/fold_0/chrombpnet_wo_bias.h5
#/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/ENCSR109RIQ/chrombpnet_model_encsr146kfx_bias/chrombpnet_wo_bias.h5
celline=fold0

chrombpnet_nb=$model
chrombpnet=$model
bias=$model
gpu=MIG-166d7783-762d-5f61-b31c-549eb4e0fba0
output_dir=$out/$celline/
interpret_dir=$out/$celline/

mkdir $output_dir
bash make_bigwig_new.sh $chrombpnet_nb $chrombpnet $bias $celline $gpu $regions1 $output_dir 



