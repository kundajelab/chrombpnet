from kerasAC.custom_losses import *
from kerasAC.metrics import *
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
                "MultichannelMultinomialNLL":MultichannelMultinomialNLL}        
from tensorflow.keras.models import load_model
model=load_model("gm12878.dnase.with.bpnet.tobias.bias.0.hdf5",custom_objects=custom_objects)
model.summary()
import pdb 
pdb.set_trace()
