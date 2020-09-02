import pdb 
from keras.models import load_model
from kerasAC.metrics import recall, specificity, fpr, fnr, precision, f1
from kerasAC.custom_losses import *
from keras.utils.generic_utils import get_custom_objects
get_custom_objects().update({"MultichannelMultinomialNLL": MultichannelMultinomialNLL})
model_string="K562.profile.peaks.only.bpnet.0.hdf5"
model=load_model(model_string)
print(model.summary())
