import pdb 
import numpy as np ;
from sklearn.metrics import average_precision_score
from kerasAC.metrics import * 
from kerasAC.custom_losses import *

import keras;

#import the various keras layers 
from keras.layers import Dense,Activation,Dropout,Flatten,Reshape,Input, Concatenate, Cropping1D
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import GlobalMaxPooling1D,MaxPooling1D
from keras.layers.normalization import BatchNormalization

from keras.optimizers import Adam
from keras.constraints import maxnorm;
from keras.regularizers import l1, l2    

from keras.models import Model

def BatchNormalization_mod(conv, bn_flag=True):
    from keras.layers.normalization import BatchNormalization
    if bn_flag:
        return BatchNormalization()(conv)
    else:
        return conv

def res_block(conv,num_filter,f_width,act,d_rate,i,bn_true=True):
    crop_id=Cropping1D(d_rate*(f_width-1))(conv)
    conv1 = BatchNormalization_mod(conv,bn_true)
    conv1 = Activation("relu")(conv1)
    conv1 = Conv1D(num_filter,f_width,dilation_rate=d_rate,padding="valid",name='conv_'+str(i)+'_a')(conv1)
    conv1 = BatchNormalization_mod(conv1,bn_true)    
    conv1 = Activation("relu")(conv1)
    conv1 = Conv1D(num_filter,f_width,dilation_rate=d_rate,padding="valid",name='conv_'+str(i)+'_b')(conv1)
    return keras.layers.Add()([conv1, crop_id])

def build1d_model_residual(input_width,input_dimension,number_of_convolutions,filters,filter_dim,dilation,activations,bn_true=True,max_flag=True):
    input1=Input(shape=(input_width,4), name='sequence')
    conv=Conv1D(32,1, padding='same',activation='relu',name = 'upsampling')(input1)
    for i in range(0,number_of_convolutions):
            conv = res_block(conv,filters[i],filter_dim[i],activations,dilation[i],i,bn_true)
    conv= Conv1D(32, 1,padding='valid', activation='relu',name='down_sampling')(conv)
    output=Conv1D(1,1,activation='relu',name='dnase')(conv)
    model = Model(inputs=[input1],outputs=[output])
    return model

def getModelGivenModelOptionsAndWeightInits():
    #read in arguments
    seed=1234
    sequence_flank=6500
    np.random.seed(seed)
    input_width=2*sequence_flank
    input_dimension=4
    number_of_convolutions=16
    filters=[32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32]
    filter_dim=[11,11,11,11,11,11,11,11,21,21,21,21,41,41,41,41]
    dilation=[1,1,1,1,4,4,4,4,10,10,10,10,25,25,25,25]
    activations='relu'
    bn_true=True

    model=build1d_model_residual(input_width,
                                 input_dimension,
                                 number_of_convolutions,
                                 filters,
                                 filter_dim,
                                 dilation,
                                 activations,
                                 bn_true=True,
                                 max_flag=True)
    
    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    loss=ambig_mean_squared_error
    model.compile(optimizer=adam,loss=loss)
    return model
