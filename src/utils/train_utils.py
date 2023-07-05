import numpy as np
import os
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping

def get_callbacks(early_stop, output_prefix):
    early_stop_callback = EarlyStopping(monitor='val_loss',
                                        patience=early_stop,
                                        verbose=1,
                                        mode='min')

    filepath = "{}.h5".format(output_prefix)
    checkpoint_callback = ModelCheckpoint(filepath,
                                          monitor='val_loss',
                                          verbose=1,
                                          save_best_only=True,
                                          mode='min')

    return [early_stop_callback, checkpoint_callback]

def get_counts_stat(cts, outputlen):
    """
    cts: N x L'
    outputlen: int L <= L'

    Compute mean reads per example in central outputlen
    """
    slice_start = cts.shape[1]//2 - outputlen//2
    slice_end = cts.shape[1]//2 + outputlen//2

    return np.mean(cts[:, slice_start:slice_end].sum(-1))
