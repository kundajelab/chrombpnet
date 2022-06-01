import os
import sys
import chrombpnet.training.utils.losses as losses
from chrombpnet.training.utils.data_utils import get_seq
import tensorflow as tf

def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"multinomial_nll": losses.mutlinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5)
    print("got the model")
    model.summary()
    return model
