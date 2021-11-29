from keras.models import load_model
import argparse
from keras.models import Model
from losses import MultichannelMultinomialNLL
import tensorflow as tf


def fetch_args():
    parser = argparse.ArgumentParser(description="remove bias component of the model")
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to trained chrombpnet model")
    parser.add_argument("-pl", "--profile_layer_name", type=str, required=True, help="Profile head layer name to use for the sequence head")
    parser.add_argument("-cl", "--counts_layer_name", type=str, required=True, help="Counts head layer name to use for the sequence head")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output model model_output_path_h5_name")
    args = parser.parse_args()
    return args

def load_model(model_hdf5):
    from keras.models import load_model
    from keras.utils.generic_utils import get_custom_objects
    custom_objects={"MultichannelMultinomialNLL": MultichannelMultinomialNLL, "tf": tf}    
    get_custom_objects().update(custom_objects)
    model=load_model(model_hdf5)
    return model

if __name__=="__main__":

    args = fetch_args()
    model=load_model(args.model_h5)
    profile_output_without_bias = model.get_layer(args.profile_layer_name).output
    counts_output_without_bias = model.get_layer(args.counts_layer_name).output
    model_without_bias = Model(inputs=model.inputs,outputs=[profile_output_without_bias, counts_output_without_bias])
    print('save model without bias') 
    model_without_bias.save(args.output_prefix+".h5")
