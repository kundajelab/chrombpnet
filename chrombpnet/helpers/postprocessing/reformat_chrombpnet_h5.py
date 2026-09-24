import argparse
import os
from keras.models import Model
from keras.layers import Input, Add, Concatenate
from chrombpnet.training.utils.layers import LogSumExp
from chrombpnet.training.utils.model_io import load_model

def parse_args():
        parser = argparse.ArgumentParser(description="Reformat chrombpnet h5 file")
        parser.add_argument("-cnb", "--chrombpnet_nb", type=str, required=True, help="Path to chrombpnet no bias model")
        parser.add_argument("-bm", "--bias_model_scaled", type=str, required=True, help="Path to scaled bias model")        
        parser.add_argument("-o", "--output_dir", type=str, required=True, help="Path to output dir")        
        args = parser.parse_args()
        return args


def chrombpnet_model(bias_model, bpnet_model_wo_bias):
	inp = Input(shape=bpnet_model_wo_bias.input_shape[1:],name='sequence')    
	bias_output=bias_model(inp)
	bpnet_model_wo_bias_new=Model(inputs=bpnet_model_wo_bias.inputs,outputs=bpnet_model_wo_bias.outputs, name="model_wo_bias")
	output_wo_bias=bpnet_model_wo_bias_new(inp)

	profile_out = Add(name="logits_profile_predictions")([output_wo_bias[0],bias_output[0]])
	concat_counts = Concatenate(axis=-1, name="concatenate")([output_wo_bias[1], bias_output[1]])
	count_out = LogSumExp(name="logcount_predictions")(concat_counts)
	model=Model(inputs=inp,outputs=[profile_out, count_out])
	return model
    
def main(args_chrombpnet_nb, args_bias, args_output_dir):

	chrombpnet_nb=load_model(args_chrombpnet_nb,compile=False)
	bias_model=load_model(args_bias,compile=False)
		
	newp = os.path.join(args_output_dir, "chrombpnet_recompiled.h5")
	new_chrom = chrombpnet_model(bias_model, chrombpnet_nb)
	new_chrom.save(newp)
	# chrombpnet 1.x also exported a TensorFlow SavedModel (chrombpnet_recompiled/); TensorFlow is no longer a
	# dependency, so only the .h5 file is written
	return newp

if __name__ == '__main__':
		
		args = parse_args()
		main(args.chrombpnet_nb, args.bias_model_scaled, args.output_dir)



