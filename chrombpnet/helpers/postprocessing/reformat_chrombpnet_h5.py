from tensorflow.keras.models import load_model
from tensorflow.keras.utils import get_custom_objects
import argparse
import tensorflow as tf
import numpy as np
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Add, Concatenate, Lambda, Flatten
import time
import os
import argparse

def parse_args():
        parser = argparse.ArgumentParser(description="Reformat chrombpnet h5 file")
        parser.add_argument("-cnb", "--chrombpnet_nb", type=str, required=True, help="Path to chrombpnet no bias model")
        parser.add_argument("-bm", "--bias_model_scaled", type=str, required=True, help="Path to scaled bias model")        
        parser.add_argument("-o", "--output_dir", type=str, required=True, help="Path to output dir")        
        args = parser.parse_args()
        return args


args = parse_args()

args_chrombpnet_nb=args.chrombpnet_nb
args_bias=args.bias_model_scaled
args_output_dir=args.output_dir

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
    
def main(args_chrombpnet_nb, args_bias, args_output_dir):

	custom_objects={"tf":tf}  
	get_custom_objects().update(custom_objects)

	chrombpnet_nb=load_model(args_chrombpnet_nb,compile=False)
	bias_model=load_model(args_bias,compile=False)
		
	newp = args_output_dir+"/chrombpnet_recompiled.h5"
	new_chrom = chrombpnet_model(bias_model, chrombpnet_nb)
	new_chrom.save(newp)
	newp = args_output_dir+"/chrombpnet_recompiled"
	new_chrom.save(newp)

if __name__ == '__main__':
		
		main(args_chrombpnet_nb, args_bias, args_output_dir)
		

	

