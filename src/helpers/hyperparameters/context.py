import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../training/')))
import utils.losses as losses
import utils.one_hot as one_hot

def load_model_wrapper(model_h5):
    # read .h5 model
    custom_objects={"MultichannelMultinomialNLL": losses.MultichannelMultinomialNLL, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_h5)
    print("got the model")
    model.summary()
    return model

