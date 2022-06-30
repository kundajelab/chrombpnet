from tensorflow import keras
from chrombpnet.training.utils import augment
from chrombpnet.training.utils import data_utils
import tensorflow as tf
import numpy as np
import random
import string
import math
import os
import json

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
    def __init__(self, peak_regions, nonpeak_regions, genome_fasta, batch_size, inputlen, outputlen, max_jitter, negative_sampling_ratio, cts_bw_file, add_revcomp, return_coords, shuffle_at_epoch_start):
        """
        seqs: B x L' x 4
        cts: B x M'
        inputlen: int (L <= L'), L' is greater to allow for cropping (= jittering)
        outputlen: int (M <= M'), M' is greater to allow for cropping (= jittering)
        batch_size: int (B)
        """

        peak_seqs, peak_cts, peak_coords, nonpeak_seqs, nonpeak_cts, nonpeak_coords, = data_utils.load_data(peak_regions, nonpeak_regions, genome_fasta, cts_bw_file, inputlen, outputlen, max_jitter)
        self.peak_seqs, self.nonpeak_seqs = peak_seqs, nonpeak_seqs
        self.peak_cts, self.nonpeak_cts = peak_cts, nonpeak_cts
        self.peak_coords, self.nonpeak_coords = peak_coords, nonpeak_coords

        self.negative_sampling_ratio = negative_sampling_ratio
        self.inputlen = inputlen
        self.outputlen = outputlen
        self.batch_size = batch_size
        self.add_revcomp = add_revcomp
        self.return_coords = return_coords
        self.shuffle_at_epoch_start = shuffle_at_epoch_start


        # random crop training data to the desired sizes, revcomp augmentation
        self.crop_revcomp_data()

    def __len__(self):

        return math.ceil(self.seqs.shape[0]/self.batch_size)


    def crop_revcomp_data(self):
        # random crop training data to inputlen and outputlen (with corresponding offsets), revcomp augmentation
        # shuffle required since otherwise peaks and nonpeaks will be together
        #Sample a fraction of the negative samples according to the specified ratio
        if (self.peak_seqs is not None) and (self.nonpeak_seqs is not None):
            # crop peak data before stacking
            cropped_peaks, cropped_cnts, cropped_coords = augment.random_crop(self.peak_seqs, self.peak_cts, self.inputlen, self.outputlen, self.peak_coords)
            #print(cropped_peaks.shape)
            #print(self.nonpeak_seqs.shape)
            if self.negative_sampling_ratio < 1.0:
                self.sampled_nonpeak_seqs, self.sampled_nonpeak_cts, self.sampled_nonpeak_coords = subsample_nonpeak_data(self.nonpeak_seqs, self.nonpeak_cts, self.nonpeak_coords, len(self.peak_seqs), self.negative_sampling_ratio)
                self.seqs = np.vstack([cropped_peaks, self.sampled_nonpeak_seqs])
                self.cts = np.vstack([cropped_cnts, self.sampled_nonpeak_cts])
                self.coords = np.vstack([cropped_coords, self.sampled_nonpeak_coords])
            else:
                self.seqs = np.vstack([cropped_peaks, self.nonpeak_seqs])
                self.cts = np.vstack([cropped_cnts, self.nonpeak_cts])
                self.coords = np.vstack([cropped_coords, self.nonpeak_coords])

        elif self.peak_seqs is not None:
            # crop peak data before stacking
            cropped_peaks, cropped_cnts, cropped_coords = augment.random_crop(self.peak_seqs, self.peak_cts, self.inputlen, self.outputlen, self.peak_coords)

            self.seqs = cropped_peaks
            self.cts = cropped_cnts
            self.coords = cropped_coords

        elif self.nonpeak_seqs is not None:
            #print(self.nonpeak_seqs.shape)

            self.seqs = self.nonpeak_seqs
            self.cts = self.nonpeak_cts
            self.coords = self.nonpeak_coords
        else :
            print("Both peak and non-peak arrays are empty")

        self.cur_seqs, self.cur_cts, self.cur_coords = augment.crop_revcomp_augment(
                                            self.seqs, self.cts, self.coords, self.inputlen, self.outputlen, 
                                            self.add_revcomp, shuffle=self.shuffle_at_epoch_start
                                          )

    def __getitem__(self, idx):
        batch_seq = self.cur_seqs[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_cts = self.cur_cts[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_coords = self.cur_coords[idx*self.batch_size:(idx+1)*self.batch_size]

        if self.return_coords:
            return (batch_seq, [batch_cts, np.log(1+batch_cts.sum(-1, keepdims=True))], batch_coords)
        else:
            return (batch_seq, [batch_cts, np.log(1+batch_cts.sum(-1, keepdims=True))])

    def on_epoch_end(self):
        self.crop_revcomp_data()

