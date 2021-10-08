from tensorflow.keras.models import load_model, model_from_json
from keras.models import Model
import tensorflow as tf

from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import *
from kerasAC.custom_losses import *



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
                        "MultichannelMultinomialNLL":MultichannelMultinomialNLL}
get_custom_objects().update(custom_objects)


#from archs import *
#model=model_from_json(open("model.0.arch","r").read())
#model.load_weights("model.0.weights")
model=load_model("/srv/scratch/anusri/chrombpnet_paper/results/chrombpnet/ATAC/K562/4_4_shifted_ATAC_09.16.2021_subsample_100M/with_k562_bias_final_model/unplug/model.0.hdf5")
print(model.summary())

#model.save('model.0.hdf5')
