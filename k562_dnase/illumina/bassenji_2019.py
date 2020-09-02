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
    #crop=Cropping1D(d_rate*(f_width-1))(conv)
    conv1 = Activation("relu")(conv)
    conv1 = Conv1D(num_filter,f_width,dilation_rate=d_rate,padding="same",name='conv_'+str(i)+'_a')(conv1)
    conv1 = BatchNormalization_mod(conv1,bn_true)
    conv1 = Activation("relu")(conv1)
    conv1 = Conv1D(768,1,padding="same",name='conv_'+str(i)+'_b')(conv1)
    conv1 = BatchNormalization_mod(conv1,bn_true)
    conv1=Dropout(0.3)(conv1)
    print('done!')
    return keras.layers.Add()([conv1, conv])

def build1d_model_residual(input_width,input_dimension,number_of_convolutions,filters,filter_dim,dilation,activations,bn_true=True,max_flag=True):
    input1=Input(shape=(input_width,4), name='sequence')
    x = Conv1D(288,15,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)
    
    x = Conv1D(339,5,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)
    
    x = Conv1D(399,5,padding='same')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(470,5,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(554,5,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(652,5,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(pool_size=2)(x)
 
    x = Conv1D(768,5,padding='valid')(input1)
    x = BatchNormalization(axis=-1)(x)
    x = Activation('relu')(x)
    conv = MaxPooling1D(pool_size=2)(x)
 
    for i in range(0,number_of_convolutions):
            conv = res_block(conv,filters[i],filter_dim[i],activations,dilation[i],i,bn_true)
    conv= Conv1D(1536, 1,padding='valid', activation='relu',name='down_sampling')(conv)
    conv=Dropout(0.05)(conv)
    output=Conv1D(1,1,activation='relu',name='dnase')(conv)
    model = Model(inputs=[input1],outputs=[output])
    return model

def getModelGivenModelOptionsAndWeightInits(args):
    #read in arguments
    seed=1234
    sequence_flank=65536
    np.random.seed(seed)
    input_width=2*sequence_flank
    input_dimension=4
    number_of_convolutions=11
    filters=[384]*number_of_convolutions
    filter_dim=[3]*number_of_convolutions
    dilation=[1,2,3,5,8,12,18,28,42,64,96]
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
    pdb.set_trace() 
    return model
