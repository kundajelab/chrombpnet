from tensorflow import keras
import archs
from utils import data_utils, train_utils, augment, argmanager
from utils.loss import multinomial_nll
import numpy as np
import random
import string
import math
import os
import json


def subsample_nonpeak_data(nonpeak_seqs, nonpeak_cts, peak_data_size, negative_sampling_ratio):
    #Randomly samples a portion of the non-peak data to use in training
    num_nonpeak_samples = int(negative_sampling_ratio * peak_data_size)
    nonpeak_indices_to_keep = np.random.choice(len(nonpeak_seqs), size=num_nonpeak_samples, replace=False)
    nonpeak_seqs = nonpeak_seqs[nonpeak_indices_to_keep]
    nonpeak_cts = nonpeak_cts[nonpeak_indices_to_keep]
    return nonpeak_seqs, nonpeak_cts


class ChromBPNetBatchGenerator(keras.utils.Sequence):
    """
    This generator randomly crops (=jitter) and revcomps training examples for 
    every epoch, and calls bias model on it, whose outputs (bias profile logits 
    and bias logcounts) are fed as input to the chrombpnet model.
    """
    def __init__(self, model_bias, peak_seqs, nonpeak_seqs, peak_cts, nonpeak_cts, negative_sampling, negative_sampling_ratio, inputlen, outputlen, batch_size):
        """
        model_bias: model that predicts bias logits and logcounts, input size L and output size M
        seqs: B x L' x 4
        cts: B x M'
        inputlen: int (L <= L'), L' is greater to allow for cropping (= jittering)
        outputlen: int (M <= M'), M' is greater to allow for cropping (= jittering)
        batch_size: int (B)
        """

        self.model_bias = model_bias
        self.peak_seqs, self.nonpeak_seqs = peak_seqs, nonpeak_seqs
        self.peak_cts, self.nonpeak_cts = peak_cts, nonpeak_cts
        self.negative_sampling = negative_sampling
        self.negative_sampling_ratio = negative_sampling_ratio
        self.inputlen = inputlen
        self.outputlen = outputlen
        self.batch_size = batch_size


        # random crop training data to the desired sizes, revcomp augmentation
        self.crop_revcomp_data()

    def __len__(self):
        return math.ceil(self.seqs.shape[0]/self.batch_size)


    def crop_revcomp_data(self):
        # random crop training data to inputlen and outputlen (with corresponding offsets), revcomp augmentation
        # shuffle required since otherwise peaks and nonpeaks will be together

        #Sample a fraction of the negative samples according to the specified ratio
        if self.negative_sampling:
            self.sampled_nonpeak_seqs, self.sampled_nonpeak_cts = subsample_nonpeak_data(self.nonpeak_seqs, self.nonpeak_cts, len(self.peak_seqs), self.negative_sampling_ratio)
            self.seqs = np.vstack([self.peak_seqs, self.sampled_nonpeak_seqs])
            self.cts = np.vstack([self.peak_cts, self.sampled_nonpeak_cts])
        else:
            self.seqs = np.vstack([self.peak_seqs, self.nonpeak_seqs])
            self.cts = np.vstack([self.peak_cts, self.nonpeak_cts])

        self.cur_seqs, self.cur_cts = augment.crop_revcomp_augment(
                                         self.seqs, self.cts, self.inputlen, self.outputlen, 
                                         shuffle=True
                                      )


        # apply bias model to those     
        self.cur_bias_logits, self.cur_bias_logcts = self.model_bias.predict(self.cur_seqs)

    def __getitem__(self, idx):
        batch_seq = self.cur_seqs[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_cts = self.cur_cts[idx*self.batch_size:(idx+1)*self.batch_size]
        
        batch_bias_logits = self.cur_bias_logits[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_bias_logcts = self.cur_bias_logcts[idx*self.batch_size:(idx+1)*self.batch_size]

        return [batch_seq, batch_bias_logits, batch_bias_logcts], [batch_cts, np.log(1+batch_cts.sum(-1, keepdims=True))] 

    def on_epoch_end(self):
        self.crop_revcomp_data()


def train_loop(model_chrombpnet, model_bias, inputlen, outputlen, train_peak_seqs, train_nonpeak_seqs, train_peak_cts, train_nonpeak_cts, 
               val_peak_seqs, val_nonpeak_seqs, val_peak_cts, val_nonpeak_cts, negative_sampling, negative_sampling_ratio, batch_size, epochs, early_stop, output_prefix): 

    if negative_sampling:
        np.random.seed(1248)
        val_nonpeak_seqs, val_nonpeak_cts = subsample_nonpeak_data(val_nonpeak_seqs, val_nonpeak_cts, len(val_peak_seqs), negative_sampling_ratio)
    val_seqs = np.vstack([val_peak_seqs, val_nonpeak_seqs])
    val_cts = np.vstack([val_peak_cts, val_nonpeak_cts])


    # need generator to crop and revcomp aug training examples, but not for 
    # validation. Also applies bias model to cropped, rev comp-ed seqs.
    train_generator = ChromBPNetBatchGenerator(model_bias, train_peak_seqs, train_nonpeak_seqs, 
                                               train_peak_cts, train_nonpeak_cts, negative_sampling, negative_sampling_ratio, inputlen, outputlen, batch_size)

    print("Train dataset length: ", len(train_generator))
    # predict validation bias logits and logcounts
    val_bias_logits, val_bias_logcts = model_bias.predict(val_seqs)

    callbacks = train_utils.get_callbacks(early_stop, output_prefix)

    history = model_chrombpnet.fit(train_generator, 
                        epochs=epochs,
                        validation_data=([val_seqs, val_bias_logits, val_bias_logcts],
                                         [val_cts, 
                                          np.log(1+val_cts.sum(-1, keepdims=True))]),
                        callbacks=callbacks)

    return history


def adjust_bias_model_logcounts(bias_model, seqs, cts):
    """
    Given a bias model, sequences and associated counts, the function adds a 
    constant to the output of the bias_model's logcounts that minimises squared
    error between predicted logcounts and observed logcounts (infered from 
    cts). This simply reduces to adding the average difference between observed 
    and predicted to the "bias" (constant additive term) of the Dense layer.

    Typically the seqs and counts would correspond to training nonpeak regions.

    ASSUMES model_bias's last layer is a dense layer that outputs logcounts. 
    This would change if you change the model.
    """

    # safeguards to prevent misuse
    assert(bias_model.layers[-1].name == "logcounts")
    assert(bias_model.layers[-1].output_shape==(None,1))
    assert(isinstance(bias_model.layers[-1], keras.layers.Dense))

    print("Predicting within adjust counts")
    _, pred_logcts = bias_model.predict(seqs, verbose=True)

    delta = np.mean(np.log(1+cts.sum(-1)) - pred_logcts.ravel())

    dw, db = bias_model.layers[-1].get_weights()
    bias_model.layers[-1].set_weights([dw, db+delta])

    return bias_model


def main():
    args = argmanager.fetch_train_chrombpnet_args()
    print(args)

    if os.path.exists("{}.h5".format(args.output_prefix)):
        raise OSError('File {}.h5 already exists'.format(args.output_prefix))

    # load bias model
    with keras.utils.CustomObjectScope({'multinomial_nll':multinomial_nll}):
        model_bias = keras.models.load_model(args.bias_model)

    # input and output shapes to bias model should be equal to inputlen and outputlen,
    # which are input and output lengths for chrombpnet model
    assert(model_bias.input_shape[1]==args.inputlen)
    assert(model_bias.output_shape[0][1]==args.outputlen)

    # load data
    train_peaks_seqs, train_peaks_cts, train_nonpeaks_seqs, train_nonpeaks_cts,\
    val_peaks_seqs, val_peaks_cts, val_nonpeaks_seqs, val_nonpeaks_cts =  \
                            data_utils.load_train_val_data(
                                args.peaks, args.nonpeaks, args.genome, args.bigwig,
                                args.val_chr, args.test_chr, args.inputlen, args.outputlen, args.max_jitter,
                                outlier=0.9999
                            )

    # adjust bias model predicted counts to accounts for sequencing depth of 
    # current sample. The adjustment is done on nonpeak regions. Note that this
    # will slightly update even if the bias model is trained on the current 
    # sample, since during training bias model, the nonpeaks regions used for
    # training are those below a stringent threshold
    if args.adjust_bias_cts:
        # crop inputlen/outputlen from center of nonpeaks
        seq_slice_start = train_nonpeaks_seqs.shape[1]//2 - args.inputlen//2
        seq_slice_end = train_nonpeaks_seqs.shape[1]//2 + args.inputlen//2
        cts_slice_start = train_nonpeaks_cts.shape[1]//2 - args.outputlen//2
        cts_slice_end = train_nonpeaks_cts.shape[1]//2 + args.outputlen//2       
        model_bias = adjust_bias_model_logcounts(model_bias, 
                                    train_nonpeaks_seqs[:, seq_slice_start:seq_slice_end], 
                                    train_nonpeaks_cts[:, cts_slice_start:cts_slice_end])

        # save adjusted bias model
        model_bias.save("{}.adjusted_bias_model.h5".format(args.output_prefix))

    # compute loss weight factor for counts loss
    counts_loss_weight = train_utils.get_counts_stat(train_peaks_cts,
                                     args.outputlen) * args.counts_weight
    print("\nCounts loss weight : {:.2f}\n".format(counts_loss_weight))

    # prepare chrombpnet model
    model_chrombpnet = archs.chrombpnet(args.inputlen, args.outputlen, args.filters, args.ndil) 
    opt = keras.optimizers.Adam(learning_rate=args.learning_rate)
    model_chrombpnet.compile(
            optimizer=opt,
            loss=[multinomial_nll, 'mse'],
            loss_weights = [1, counts_loss_weight]
        )

    history = train_loop(model_chrombpnet, model_bias, args.inputlen, args.outputlen, 
                         train_peaks_seqs, train_nonpeaks_seqs,
                         train_peaks_cts, train_nonpeaks_cts,
                         val_peaks_seqs, val_nonpeaks_seqs,
                         val_peaks_cts, val_nonpeaks_cts, args.negative_sampling, args.negative_sampling_ratio,
                         args.batch_size, args.epochs, 
                         args.early_stop, args.output_prefix)

    with open("{}.history.json".format(args.output_prefix), "w") as f:
        json.dump(history.history, f, ensure_ascii=False, indent=4)

if __name__=="__main__":
    main()

