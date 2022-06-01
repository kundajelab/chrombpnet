# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

python tf_modiscohits.py --outdir=count/k562_hits/ \
	/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_nov_08/ATAC/K562/4_4_shifted_ATAC_09.29.2021_bias_filters_500/final_model_step3/unplug/deepshap/20K.fold0.deepSHAP \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/SIGNAL/K562/4_4_shifted_ATAC_09.29.2021_bias_filters_500_new/20K.fold0.deepSHAP.count_shap.hdf5 \
	/srv/scratch/anusri/bassik_controls/output_summits/final.bed \


python tf_modiscohits.py --outdir=count/gm12878_hits/ \
	/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_nov_08/ATAC/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500/final_model_step3/unplug/deepshap/20K.fold0.deepSHAP \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/SIGNAL/GM12878/4_4_shifted_ATAC_09.06.2021_bias_filters_500_new/20K.fold0.deepSHAP.count_shap.hdf5 \
	/srv/scratch/anusri/bassik_controls/output_summits/final.bed \


python tf_modiscohits.py --outdir=count/hepg2_hits/ \
	/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_sept_09/HEPG2/ATAC_07.22.2021/with_gm12878_bias_final_model/unplug/deepshap/20K.fold0.deepSHAP \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/SIGNAL/HEPG2/ATAC_07.22.2021_withgm12878bias/20K.fold0.deepSHAP.count_shap.hdf5 \
	/srv/scratch/anusri/bassik_controls/output_summits/final.bed \

python tf_modiscohits.py --outdir=count/imr90_hits/ \
	/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_sept_09/IMR90/ATAC_07.22.2021/with_gm12878_bias_final_model/unplug/deepshap/20K.fold0.deepSHAP \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/SIGNAL/IMR90/ATAC_07.22.2021_withgm12878bias/20K.fold0.deepSHAP.count_shap.hdf5 \
	/srv/scratch/anusri/bassik_controls/output_summits/final.bed \

python tf_modiscohits.py --outdir=count/h1_hits/ \
	/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet_sept_09/H1/ATAC_07.22.2021/with_gm12878_bias_final_model/unplug/deepshap/20K.fold0.deepSHAP \
	/oak/stanford/groups/akundaje/projects/chrombpnet_paper/modisco/flank_len_500/SIGNAL/H1/ATAC_07.22.2021_withgm12878bias/20K.fold0.deepSHAP.count_shap.hdf5 \
	/srv/scratch/anusri/bassik_controls/output_summits/final.bed \



