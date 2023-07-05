import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from . import one_hot

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def get_cts(peaks_df, bw, width):
    """
    Fetches values from a bigwig bw, given a df with minimally
    chr, start and summit columns. Summit is relative to start.
    Retrieves values of specified width centered at summit.

    "cts" = per base counts across a region
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append(np.nan_to_num(bw.values(r['chr'], 
                                            r['start'] + r['summit'] - width//2,
                                            r['start'] + r['summit'] + width//2)))
        
    return np.array(vals)

def get_seq(peaks_df, genome, width):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append(str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)]))

    return one_hot.dna_to_one_hot(vals)

def get_seq_cts(peaks_df, genome, bw, input_width, output_width):
    seq = get_seq(peaks_df, genome, input_width)
    cts = get_cts(peaks_df, bw, output_width)

    return seq, cts

def load_train_val_data(peaks_file, nonpeaks_file, genome_fasta, cts_bw_file, 
                        val_chrs, test_chrs, inputlen, outputlen, max_jitter, outlier=0.9999):
    """
    Load sequences and corresponding base resolution counts for training, 
    validation regions in peaks and nonpeaks (2 x 2 x 2 = 8 matrices).

    For training peaks/nonpeaks, values for inputlen + 2*max_jitter and outputlen + 2*max_jitter 
    are returned centered at peak summit. This allows for jittering examples by randomly
    cropping. Data of width inputlen/outputlen is returned for validation
    data.

    If outliers is not None, removes training examples with counts > outlier%ile
    """
    assert(inputlen%2==0)
    assert(outputlen%2==0)

    peaks_df = pd.read_csv(peaks_file, sep='\t', names=NARROWPEAK_SCHEMA)
    nonpeaks_df = pd.read_csv(nonpeaks_file, sep='\t', names=NARROWPEAK_SCHEMA)

    train_peaks_df = peaks_df[~peaks_df['chr'].isin(test_chrs+val_chrs)]
    val_peaks_df = peaks_df[peaks_df['chr'].isin(val_chrs)]

    train_nonpeaks_df = nonpeaks_df[~nonpeaks_df['chr'].isin(test_chrs+val_chrs)]
    val_nonpeaks_df = nonpeaks_df[nonpeaks_df['chr'].isin(val_chrs)]

    cts_bw = pyBigWig.open(cts_bw_file)
    genome = pyfaidx.Fasta(genome_fasta)

    train_peaks_seqs, train_peaks_cts = get_seq_cts(train_peaks_df,
                                                   genome,
                                                   cts_bw,
                                                   inputlen+2*max_jitter,
                                                   outputlen+2*max_jitter)

    train_nonpeaks_seqs, train_nonpeaks_cts = get_seq_cts(train_nonpeaks_df,
                                                          genome,
                                                          cts_bw,
                                                          inputlen+2*max_jitter,
                                                          outputlen+2*max_jitter) 
    
    val_peaks_seqs, val_peaks_cts = get_seq_cts(val_peaks_df,
                                               genome,
                                               cts_bw,
                                               inputlen,
                                               outputlen)

    val_nonpeaks_seqs, val_nonpeaks_cts = get_seq_cts(val_nonpeaks_df,
                                                      genome,
                                                      cts_bw,
                                                      inputlen,
                                                      outputlen)

    cts_bw.close()
    genome.close()

    if outlier is not None:
        thresh = np.quantile(np.vstack([train_peaks_cts, train_nonpeaks_cts]).sum(-1), outlier)
        train_peaks_subset = train_peaks_cts.sum(-1) < thresh
        train_nonpeaks_subset = train_nonpeaks_cts.sum(-1) < thresh
       
        train_peaks_seqs = train_peaks_seqs[train_peaks_subset]
        train_peaks_cts = train_peaks_cts[train_peaks_subset]
        train_nonpeaks_seqs = train_nonpeaks_seqs[train_nonpeaks_subset]
        train_nonpeaks_cts = train_nonpeaks_cts[train_nonpeaks_subset]

    return(train_peaks_seqs, train_peaks_cts, train_nonpeaks_seqs, train_nonpeaks_cts, 
           val_peaks_seqs, val_peaks_cts, val_nonpeaks_seqs, val_nonpeaks_cts)


def load_test_data(peaks_file, nonpeaks_file, genome_fasta, cts_bw_file,
                   test_chrs, inputlen, outputlen):
    assert(inputlen%2==0)

    peaks_df = pd.read_csv(peaks_file, sep='\t', names=NARROWPEAK_SCHEMA)
    nonpeaks_df = pd.read_csv(nonpeaks_file, sep='\t', names=NARROWPEAK_SCHEMA)

    test_peaks_df = peaks_df[peaks_df['chr'].isin(test_chrs)]
    test_nonpeaks_df = nonpeaks_df[nonpeaks_df['chr'].isin(test_chrs)]

    cts_bw = pyBigWig.open(cts_bw_file)
    genome = pyfaidx.Fasta(genome_fasta)

    test_peaks_seqs, test_peaks_cts = get_seq_cts(test_peaks_df,
                                                  genome,
                                                  cts_bw,
                                                  inputlen,
                                                  outputlen)

    test_nonpeaks_seqs, test_nonpeaks_cts = get_seq_cts(test_nonpeaks_df,
                                                        genome,
                                                        cts_bw,
                                                        inputlen,
                                                        outputlen)

    cts_bw.close()
    genome.close()

    return test_peaks_seqs, test_peaks_cts, test_nonpeaks_seqs, test_nonpeaks_cts
