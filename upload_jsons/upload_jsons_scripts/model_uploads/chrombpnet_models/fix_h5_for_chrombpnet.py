from tensorflow.keras.models import load_model
from tensorflow.keras.utils import get_custom_objects
import argparse
import tensorflow as tf
import numpy as np
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Add, Concatenate, Lambda, Flatten
import time
import os

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"
primary_models_path = ["chrombppnet_model_encsr283tme_bias", "chrombppnet_model_encsr283tme_bias_fold_1", "chrombppnet_model_encsr283tme_bias_fold_2", "chrombppnet_model_encsr283tme_bias_fold_3", "chrombppnet_model_encsr283tme_bias_fold_4"]
celline_models_path = ["chrombpnet_model_feb15_fold_0", "chrombpnet_model_feb15_fold_1", "chrombpnet_model_feb15_fold_2", "chrombpnet_model_feb15_fold_3", "chrombpnet_model_feb15_fold_4"]
tissue_models_path = ["chrombpnet_model_encsr880cub_bias","chrombppnet_model_encsr880cub_bias_fold_1","chrombppnet_model_encsr880cub_bias_fold_2","chrombppnet_model_encsr880cub_bias_fold_3","chrombppnet_model_encsr880cub_bias_fold_4"]
invitro_models_path = ["chrombpnet_model_encsr146kfx_bias", "chrombpnet_model_encsr146kfx_bias_fold_1", "chrombpnet_model_encsr146kfx_bias_fold_2", "chrombpnet_model_encsr146kfx_bias_fold_3", "chrombpnet_model_encsr146kfx_bias_fold_4"]


tissue_encids = open("data/tissue_passed.txt").readlines()
tissue_encids = [line.strip() for line in tissue_encids]

primary_encids = open("data/primary_passed.txt").readlines()
primary_encids = [line.strip() for line in primary_encids]

celline_encids = open("data/cellline_passed.txt").readlines()
celline_encids = [line.strip() for line in celline_encids]

invitro_encids = open("data/invitro_passed.txt").readlines()
invitro_encids = [line.strip() for line in invitro_encids]


def chrombpnet_model(bias_model, bpnet_model_wo_bias):
	inp = Input(shape=(2114, 4),name='sequence')    
	bias_output=bias_model(inp)
	bpnet_model_wo_bias_new=Model(inputs=bpnet_model_wo_bias.inputs,outputs=bpnet_model_wo_bias.outputs, name="model_wo_bias")
	output_wo_bias=bpnet_model_wo_bias_new(inp)

	profile_out = Add(name="logits_profile_predictions")([output_wo_bias[0],bias_output[0]])
	concat_counts = Concatenate(axis=-1)([output_wo_bias[1], bias_output[1]])
	count_out = Lambda(lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
						name="logcount_predictions")(concat_counts)
	model=Model(inputs=[inp],outputs=[profile_out, count_out])
	return model
    
def main(args_chrombpnet_nb, args_chrombpnet, args_bias, args_output_dir):

	custom_objects={"tf":tf}  
	get_custom_objects().update(custom_objects)
	
	if os.path.isfile(args_chrombpnet):

		chrombpnet_nb=load_model(args_chrombpnet_nb,compile=False)
		chrombpnet=load_model(args_chrombpnet,compile=False)
		bias_model=load_model(args_bias,compile=False)

		model2 = chrombpnet_nb
		temp1 = chrombpnet.get_layer("model_wo_bias").output
		model1 = Model(inputs=chrombpnet.get_layer("model_wo_bias").inputs,outputs=[temp1[0], temp1[1]])

		print(chrombpnet_nb.summary())
		layer_names=[layer.name for layer in chrombpnet_nb.layers]
		layer_names=["wo_bias_bpnet_1st_conv" ,"wo_bias_bpnet_1conv", "wo_bias_bpnet_2conv", "wo_bias_bpnet_3conv", "wo_bias_bpnet_4conv", "wo_bias_bpnet_5conv", "wo_bias_bpnet_6conv", "wo_bias_bpnet_7conv", "wo_bias_bpnet_8conv", "wo_bias_bpnet_logcount_predictions", "wo_bias_bpnet_prof_out_precrop"]
		checks = []
		for nl in layer_names:
			print(nl)
			weight1 = chrombpnet.get_layer("model_wo_bias").get_layer(nl).weights[0]
			weight2 = chrombpnet_nb.get_layer(nl).weights[0]
			checks.append(np.array_equal(weight1, weight2))

			weight1 = chrombpnet.get_layer("model_wo_bias").get_layer(nl).weights[1]
			weight2 = chrombpnet_nb.get_layer(nl).weights[1]
			checks.append(np.array_equal(weight1, weight2))

		if len(checks)==sum(checks):
			print(True)
			f = open(args_output_dir+"/check_passed.txt", "w")
			f.write("done")
			f.close()
		else:
			newp = args_output_dir+"/chrombpnet_new.h5"
			new_chrom = chrombpnet_model(bias_model, chrombpnet_nb)
			new_chrom.save(newp)
			newp = args_output_dir+"/chrombpnet_new"
			new_chrom.save(newp)
			print(False)
	else:
	
		chrombpnet_nb=load_model(args_chrombpnet_nb,compile=False)
		bias_model=load_model(args_bias,compile=False)
		
		newp = args_output_dir+"/chrombpnet_new.h5"
		new_chrom = chrombpnet_model(bias_model, chrombpnet_nb)
		new_chrom.save(newp)
		newp = args_output_dir+"/chrombpnet_new"
		new_chrom.save(newp)
		print(False)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-e", "--encid", type=str)
	args = parser.parse_args()
	encid=args.encid
	
	if encid in primary_encids:
		models_path = primary_models_path
		print("primary")
	elif encid in tissue_encids:
		models_path = tissue_models_path
		print("tissue")
	elif encid in invitro_encids:
		models_path = invitro_models_path
		print("invitro")
	elif encid in celline_encids:
		models_path = celline_models_path
		print("celline")
	else:
		print(encid)
	
	for idx in range(5):
		args_chrombpnet_nb=os.path.join(odir+encid,models_path[idx]+"/chrombpnet_wo_bias.h5")
		args_chrombpnet=os.path.join(odir+encid,models_path[idx]+"/chrombpnet.h5")
		args_bias=os.path.join(odir+encid,models_path[idx]+"/bias_model_scaled.h5")
		args_output_dir=os.path.join(odir+encid,models_path[idx]+"/new_chrombpnet_model/")
		if not os.path.isdir(args_output_dir):
			os.mkdir(args_output_dir)
		main(args_chrombpnet_nb, args_chrombpnet, args_bias, args_output_dir)
		

	

