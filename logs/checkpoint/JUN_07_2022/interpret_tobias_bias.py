import pyBigWig
import numpy as np
import matplotlib
matplotlib.use('Agg')

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

import matplotlib.pyplot as plt
import h5py
import tensorflow as tf
import src.training.utils.losses as losses
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
def load_model_wrapper(model_h5):
    # read .h5 model
    custom_objects={"tf": tf, "multinomial_nll":losses.multinomial_nll}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_h5)
    print("got the model")
    model.summary()
    return model

bw = pyBigWig.open("results/tobias/ATAC_PE/GM12878/data/GM12878_expected_not_multiplied_with_counts/sorted_merged_expected.bw")
data = bw.values('chr1', 100027296+111-500-50, 100027296+111+500+50)
datav = np.nan_to_num(data)

model="results/tobias/ATAC_PE/GM12878/ATAC_PE_03.11.2022_tobias/tobias_bias_model/tobias_bias.h5"
model=load_model_wrapper(model)
print(np.array(np.expand_dims([datav], axis=-1)).shape)
print(np.log(np.sum(datav)+1).shape)
pred_m = model.predict_on_batch([np.array(np.expand_dims([datav], axis=-1)), np.expand_dims(np.log(np.sum(datav)+1),axis=-1)])


print(softmax(pred_m[0]))
print(pred_m[1])
#plt.plot(datav/np.sum(datav))
plt.plot(softmax(pred_m[0])[0])

plt.savefig("temp.png")
