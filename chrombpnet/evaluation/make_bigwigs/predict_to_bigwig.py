import argparse
import pyBigWig
import numpy as np
import pandas as pd
import pyfaidx
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import tensorflow as tf
import chrombpnet.evaluation.make_bigwigs.bigwig_helper as bigwig_helper
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.data_utils as data_utils 
import chrombpnet.training.utils.one_hot as one_hot

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


# need full paths!
def parse_args():
    parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS!")
    parser.add_argument("-bm", "--bias-model", type=str, required=True, help="Path to bias model h5")
    parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model h5")
    parser.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=True, help="Path to chrombpnet model h5")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the input regions are assumed to be for [start+summit-(inputlen/2):start+summit+(inputlen/2)]. Should not be piped since it is read twice!")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
    parser.add_argument("-o", "--out-prefix", type=str, required=True, help="Output bigwig file")
    parser.add_argument("-b", "--batch-size", type=int, default=64)
    parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
    parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    args = parser.parse_args()
    print(args)
    return args

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def load_model_wrapper(model_hdf5):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_hdf5)
    print("got the model")
    model.summary()
    return model

def main():
    args=parse_args()
    model_chrombpnet_nb = load_model_wrapper(model_hdf5=args.chrombpnet_model_nb)
    model_chrombpnet = load_model_wrapper(model_hdf5=args.chrombpnet_model)
    model_bias = load_model_wrapper(model_hdf5=args.bias_model)


    inputlen = int(model_bias.input_shape[1])
    outputlen = int(model_bias.output_shape[0][1])
    # input and output shapes should be the same for bias model and
    # chrombpnet model
    assert(model_chrombpnet.input_shape[1]==inputlen)
    assert(model_chrombpnet.output_shape[0][1]==outputlen)

    # load data
    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
    print(regions_df.head())
    gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
    regions = bigwig_helper.get_regions(args.regions, outputlen) # output regions

    if args.debug_chr is not None:
        regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
        regions = [x for x in regions if x[0]==args.debug_chr]

    with pyfaidx.Fasta(args.genome) as g:
        seqs = data_utils.get_seq(regions_df, g, inputlen)

    pred_bias_logits, pred_bias_logcts = model_bias.predict(seqs,
                                          batch_size = args.batch_size,
                                          verbose=True)
    pred_logits, pred_logcts = model_chrombpnet.predict([seqs],
                                          batch_size = args.batch_size,
                                          verbose=True)
    pred_logits_wo_bias, pred_logcts_wo_bias = model_chrombpnet_nb.predict([seqs],
                                          batch_size = args.batch_size,
                                          verbose=True)

    pred_bias_logits = np.squeeze(pred_bias_logits)
    pred_logits = np.squeeze(pred_logits)

    pred_logits_wo_bias = np.squeeze(pred_logits_wo_bias)
    bigwig_helper.write_bigwig(softmax(pred_bias_logits) * (np.expand_dims(np.exp(pred_bias_logcts)[:,0],axis=1)), 
                               regions, 
                               gs, 
                               args.out_prefix + "bias.bw", 
                               args.out_prefix + "bias.stats.txt", 
                               debug_chr=args.debug_chr, 
                               use_tqdm=args.tqdm)

    bigwig_helper.write_bigwig(softmax(pred_logits_wo_bias) * (np.expand_dims(np.exp(pred_logcts_wo_bias)[:,0],axis=1)), 
                               regions, 
                               gs, 
                               args.out_prefix + "_wo_bias.bw", 
                               args.out_prefix + "_wo_bias.stats.txt", 
                               debug_chr=args.debug_chr, 
                               use_tqdm=args.tqdm)

    bigwig_helper.write_bigwig(softmax(pred_logits) * (np.expand_dims(np.exp(pred_logcts)[:,0],axis=1)),
                               regions,
                               gs,
                               args.out_prefix + "_w_bias.bw",
                               args.out_prefix + "_w_bias.stats.txt",
                               debug_chr=args.debug_chr,
                               use_tqdm=args.tqdm)
    
