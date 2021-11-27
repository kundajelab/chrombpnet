from tensorflow import keras
import utils.batchgen_generator_augment as augment
import utils.batchgen_generator_utils as utils
import tensorflow as tf
import numpy as np
import random
import string
import math
import os
import json

os.environ['PYTHONHASHSEED'] = '0'

def subsample_nonpeak_data(nonpeak_seqs, nonpeak_cts, nonpeak_coords, peak_data_size, negative_sampling_ratio):
    #Randomly samples a portion of the non-peak data to use in training
    num_nonpeak_samples = int(negative_sampling_ratio * peak_data_size)
    nonpeak_indices_to_keep = np.random.choice(len(nonpeak_seqs), size=num_nonpeak_samples, replace=False)
    nonpeak_seqs = nonpeak_seqs[nonpeak_indices_to_keep]
    nonpeak_cts = nonpeak_cts[nonpeak_indices_to_keep]
    nonpeak_coords = nonpeak_coords[nonpeak_indices_to_keep]
    return nonpeak_seqs, nonpeak_cts, nonpeak_coords

class ChromBPNetBatchGenerator(keras.utils.Sequence):
    """
    This generator randomly crops (=jitter) and revcomps training examples for 
    every epoch, and calls bias model on it, whose outputs (bias profile logits 
    and bias logcounts) are fed as input to the chrombpnet model.
    """
    def __init__(self, peak_regions, nonpeak_regions, genome_fasta, batch_size, inputlen, outputlen, max_jitter, negative_sampling_ratio, cts_sum_min_thresh, cts_sum_max_thresh, cts_bw_file, seed, add_revcomp, return_coords):
        """
        seqs: B x L' x 4
        cts: B x M'
        inputlen: int (L <= L'), L' is greater to allow for cropping (= jittering)
        outputlen: int (M <= M'), M' is greater to allow for cropping (= jittering)
        batch_size: int (B)
        """

        np.random.seed(seed)
        random.seed(seed)
        tf.random.set_seed(seed)

        peak_seqs, peak_cts, peak_coords, nonpeak_seqs, nonpeak_cts, nonpeak_coords, = utils.load_data(peak_regions, nonpeak_regions, genome_fasta, cts_bw_file, inputlen, outputlen, max_jitter, cts_sum_min_thresh, cts_sum_max_thresh)
        self.peak_seqs, self.nonpeak_seqs = peak_seqs, nonpeak_seqs
        self.peak_cts, self.nonpeak_cts = peak_cts, nonpeak_cts
        self.peak_coords, self.nonpeak_coords = peak_coords, nonpeak_coords

        self.negative_sampling_ratio = negative_sampling_ratio
        self.inputlen = inputlen
        self.outputlen = outputlen
        self.batch_size = batch_size
        self.add_revcomp = add_revcomp
        self.return_coords = return_coords
        self.seed=seed


        # random crop training data to the desired sizes, revcomp augmentation
        self.crop_revcomp_data()

    def __len__(self):
        return math.ceil(self.seqs.shape[0]/self.batch_size)


    def crop_revcomp_data(self):
        # random crop training data to inputlen and outputlen (with corresponding offsets), revcomp augmentation
        # shuffle required since otherwise peaks and nonpeaks will be together
        #Sample a fraction of the negative samples according to the specified ratio

        if (self.peak_seqs is not None) and (self.nonpeak_seqs is not None):
            if self.negative_sampling_ratio < 1.0:
                self.sampled_nonpeak_seqs, self.sampled_nonpeak_cts, self.sampled_nonpeak_coords = subsample_nonpeak_data(self.nonpeak_seqs, self.nonpeak_cts, self.nonpeak_coords, len(self.peak_seqs), self.negative_sampling_ratio)
                self.seqs = np.vstack([self.peak_seqs, self.sampled_nonpeak_seqs])
                self.cts = np.vstack([self.peak_cts, self.sampled_nonpeak_cts])
                self.coords = np.vstack([self.peak_coords, self.sampled_nonpeak_coords])
            else:
                self.seqs = np.vstack([self.peak_seqs, self.nonpeak_seqs])
                self.cts = np.vstack([self.peak_cts, self.nonpeak_cts])
                self.coords = np.vstack([self.peak_coords, self.nonpeak_coords])
        elif self.peak_seqs is not None:
            self.seqs = self.peak_seqs
            self.cts = self.peak_cts
            self.coords = self.peak_coords
        elif self.nonpeak_seqs is not None:
            self.seqs = self.nonpeak_seqs
            self.cts = self.nonpeak_cts
            self.coords = self.nonpeak_coords
        else :
            print("Both peak and non-peak arrays are empty")


        if self.add_revcomp:
            self.cur_seqs, self.cur_cts, self.cur_coords = augment.crop_revcomp_augment(
                                             self.seqs, self.cts, self.coords, self.inputlen, self.outputlen, 
                                             seed=self.seed, shuffle=True
                                          )
        else:
            self.cur_seqs, self.cur_cts, self.cur_coords = self.seqs, self.cts, self.coords


    def __getitem__(self, idx):
        batch_seq = self.cur_seqs[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_cts = self.cur_cts[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_coords = self.cur_coords[idx*self.batch_size:(idx+1)*self.batch_size]
        
        if self.return_coords:
            return (batch_seq, [np.expand_dims(batch_cts,axis=2), np.log(1+batch_cts.sum(-1, keepdims=True))], batch_coords)
        else:
            return (batch_seq, [np.expand_dims(batch_cts,axis=2), np.log(1+batch_cts.sum(-1, keepdims=True))])

    def on_epoch_end(self):
        self.crop_revcomp_data()

