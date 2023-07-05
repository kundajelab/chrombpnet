from tensorflow import keras
import archs
import numpy as np
from utils import data_utils, train_utils, augment, argmanager 
from utils.loss import multinomial_nll
import random
import string
import math
import json
import os


class BiasBatchGenerator(keras.utils.Sequence):
    """
    This generator randomly crops (=jitter) and revcomps training examples 
    for every epoch.
    """
    def __init__(self, seqs, cts, inputlen, outputlen, batch_size):
        """
        seqs: B x L' x 4
        cts: B x M'
        inputlen: int (L <= L'), L' is greater to allow for cropping (= jittering)
        outputlen: int (M <= M'), M' is greater to allow for cropping (= jittering)
        batch_size: int (B)
        """
        self.seqs = seqs
        self.cts = cts
        self.inputlen = inputlen
        self.outputlen = outputlen
        self.batch_size = batch_size

        # random crop training data to appropriate length, revcomp augmentation
        self.crop_revcomp_data()

    def __len__(self):
        return math.ceil(self.seqs.shape[0]/self.batch_size)

    def __getitem__(self, idx):
        batch_seq = self.cur_seqs[idx*self.batch_size:(idx+1)*self.batch_size]
        batch_cts = self.cur_cts[idx*self.batch_size:(idx+1)*self.batch_size]

        return batch_seq, [batch_cts , np.log(1+batch_cts.sum(-1, keepdims=True))]

    def crop_revcomp_data(self):
        self.cur_seqs, self.cur_cts = augment.crop_revcomp_augment(
                                         self.seqs, self.cts, self.inputlen, self.outputlen,
                                         shuffle=True
                                      )

    def on_epoch_end(self):
        # random crop training data to inputlen (seq) and outputlen (cts), revcomp augmentation
        self.crop_revcomp_data()


def train_loop(model, inputlen, outputlen, train_seqs, train_cts, val_seqs, val_cts, 
               batch_size, epochs, early_stop, output_prefix):

    # need generator to crop and revcomp aug training examples, but not for validation
    train_generator = BiasBatchGenerator(train_seqs, train_cts, inputlen, outputlen, batch_size)

    callbacks = train_utils.get_callbacks(early_stop, output_prefix)

    history = model.fit(train_generator, 
                        epochs=epochs,
                        validation_data=(val_seqs, 
                                         [val_cts, 
                                          np.log(1+val_cts.sum(-1, keepdims=True))]),
                        callbacks=callbacks)

    return history


def compute_threshold(train_peak_cts, threshold_factor):
    """
    Not all nonpeaks reigons necessarily have low counts. This method computes
    a stringent threshold based on total counts in training peak regions.
    """
    return train_peak_cts.sum(-1).min() * threshold_factor


def filter_nonpeaks_dat(train_nonpeaks_seqs, train_nonpeaks_cts, val_nonpeaks_seqs, 
                        val_nonpeaks_cts, max_totcount_thresh):
    train_select = train_nonpeaks_cts.sum(-1) <= max_totcount_thresh
    val_select = val_nonpeaks_cts.sum(-1) <= max_totcount_thresh
    
    print("\nMax counts threshold is {:d}. Retaining {:d} ({:.2f}%) training examples and {:d} ({:.2f}%) validation examples".format(
        int(max_totcount_thresh),
        int(np.sum(train_select)), 100*np.mean(train_select), 
        int(np.sum(val_select)), 100*np.mean(val_select)))

    return train_nonpeaks_seqs[train_select], train_nonpeaks_cts[train_select], val_nonpeaks_seqs[val_select], val_nonpeaks_cts[val_select]


def main():
    args = argmanager.fetch_train_bias_args()
    print(args)

    if os.path.exists("{}.h5".format(args.output_prefix)):
        raise OSError('File {}.h5 already exists'.format(args.output_prefix))

    _, train_peak_cts, train_nonpeaks_seqs, train_nonpeaks_cts,\
    _, _, val_nonpeaks_seqs, val_nonpeaks_cts =  data_utils.load_train_val_data(
                                args.peaks, args.nonpeaks, args.genome, args.bigwig,
                                args.val_chr, args.test_chr, args.inputlen, args.outputlen, args.max_jitter,
                                outlier=None)

    # filter nonpeaks to be below a certain stringent total count threshold
    max_totcount = compute_threshold(train_peak_cts, args.threshold_factor)

    train_nonpeaks_seqs, train_nonpeaks_cts, \
    val_nonpeaks_seqs, val_nonpeaks_cts  = filter_nonpeaks_dat(
                                    train_nonpeaks_seqs, train_nonpeaks_cts,
                                    val_nonpeaks_seqs, val_nonpeaks_cts,
                                    max_totcount)

    # compute loss weight factor for counts
    counts_loss_weight = train_utils.get_counts_stat(
                             train_nonpeaks_cts, 
                             args.outputlen) * args.counts_weight
    print("\nCounts loss weight : {:.2f}\n".format(counts_loss_weight))

    # prepare bias bpnet model
    model = archs.bpnet_seq(args.inputlen, args.outputlen, args.filters, args.ndil) 
    opt = keras.optimizers.Adam(learning_rate=args.learning_rate)

    model.compile(
            optimizer=opt,
            loss=[multinomial_nll, 'mse'],
            loss_weights = [1, counts_loss_weight]
        )

    history = train_loop(model, args.inputlen, args.outputlen, train_nonpeaks_seqs, 
                         train_nonpeaks_cts, val_nonpeaks_seqs, 
                         val_nonpeaks_cts, args.batch_size, args.epochs, 
                         args.early_stop, args.output_prefix)

    with open("{}.history.json".format(args.output_prefix), "w") as f:
        json.dump(history.history, f, ensure_ascii=False, indent=4)

if __name__=="__main__":
    main()

