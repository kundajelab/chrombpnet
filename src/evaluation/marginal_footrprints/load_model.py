import pdb 
import argparse
import pickle
import tensorflow
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import kerasAC 
from kerasAC.generators.tiledb_predict_generator import *
from kerasAC.tiledb_config import *
from kerasAC.interpret.deepshap import * 
from kerasAC.interpret.profile_shap import * 
from kerasAC.helpers.transform_bpnet_io import *
#load the model!
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import * 
from kerasAC.custom_losses import * 
from tensorflow.keras.models import load_model, model_from_json
#from bpnet_shap_weighted import *
from keras_genomics.layers.convolutional import RevCompConv1D
import keras
def load_model_weights(weight_file,model):
    if weight_file is None:
        #nothing to do
        return model
    #sync the model locally if it's on AWS
    if weight_file.startswith("s3://"):
        s3_model_weights=download_s3_file(weight_file)
    else:
        s3_model_weights=weight_file
    import  h5py
    try:
        model.load_weights(s3_model_weights)#,by_name=True)
    except:
        with h5py.File(s3_model_weights,'r') as file:
            weight_file=file['model_1']
            for layer in model.layers:
                try:
                    layer_weights=weight_file[layer.name]
                except:
                    print('no weights files saved for layer:'+str(layer.name))
                    continue
                try:
                    weights = []
                    # Extract weights
                    for term in layer_weights:
                        if isinstance(layer_weights[term], h5py.Dataset):
                            # Convert weights to numpy array and prepend to list
                            weights.insert(0, np.array(layer_weights[term]))
                    # Load weights to model
                    layer.set_weights(weights)
                    print("loaded weights for layer:"+str(layer.name))
                except Exception as e:
                    print("Error: Could not load weights for layer:"+str(layer.name))
                    raise e
    return model

class BiasLayer(keras.layers.Layer):
    def __init__(self, *args, **kwargs):
        super(BiasLayer, self).__init__(*args, **kwargs)

    def build(self, input_shape):
        self.bias = self.add_weight('bias',
                                    shape=input_shape[1:],
                                    initializer='zeros',
                                    trainable=True)
    def call(self, x):
        return x + self.bias

def load_model_wrapper(model_hdf5=None,json_string=None, weights=None): 
    custom_objects={"recall":recall,
                        "sensitivity":recall,
                        "specificity":specificity,
                        "fpr":fpr,
                        "fnr":fnr,
                        "precision":precision,
                        "f1":f1,
                        "ambig_binary_crossentropy":ambig_binary_crossentropy,
                        "ambig_mean_absolute_error":ambig_mean_absolute_error,
                        "ambig_mean_squared_error":ambig_mean_squared_error,
                        "MultichannelMultinomialMSE":MultichannelMultinomialMSE,
                        "RevCompConv1D":RevCompConv1D,
                        "BiasLayer": BiasLayer,
                        "MultichannelMultinomialNLL":MultichannelMultinomialNLL}
    get_custom_objects().update(custom_objects)
    if model_hdf5:
        model=load_model(model_hdf5)
        print(model.summary())
    else:
        model=model_from_json(open(json_string,"r").read())
        model=load_model_weights(weights,model)
    return model

