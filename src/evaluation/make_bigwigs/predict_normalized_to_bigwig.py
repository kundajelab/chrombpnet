import tensorflow as tfc
from tensorflow import keras
import argparse
import pyBigWig
import numpy as np
import pandas as pd
import bigwig_helper
import pyfaidx
import sys
from context import data_utils as data_utils
from context import load_model_wrapper as load_model_wrapper
from sklearn.isotonic import IsotonicRegression
import json
NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


# need full paths!
parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS! For now, this script needs to be run from within current directory. ")
parser.add_argument("-bm", "--bias-model", type=str, required=True, help="Path to bias model h5")
parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model h5")
parser.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=True, help="Path to chrombpnet model h5")
parser.add_argument("-obs", "--observed", type=str, required=True, help="Path to observed")
parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the input regions are assumed to be for [start+summit-(inputlen/2):start+summit+(inputlen/2)]. Should not be piped since it is read twice!")
parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
parser.add_argument("-o", "--out-prefix", type=str, required=True, help="Output bigwig file")
parser.add_argument("-b", "--batch-size", type=int, default=64)
parser.add_argument("-f", "--fold", type=str, required=False)
parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")

args = parser.parse_args()


print(args)

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)


model_chrombpnet_nb = load_model_wrapper(model_hdf5=args.chrombpnet_model_nb)
model_chrombpnet = load_model_wrapper(model_hdf5=args.chrombpnet_model)
model_bias = load_model_wrapper(model_hdf5=args.bias_model)

folds = json.load(open(args.fold))
valid_chroms=folds["valid"]

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


bw = pyBigWig.open(args.observed)
outs = data_utils.get_cts(regions_df[regions_df['chr'].isin(valid_chroms)], bw, outputlen)
outs_obs = data_utils.get_cts(regions_df, bw, outputlen)
print(outs.shape)
observed_log_cnts = np.log(np.sum(outs,axis=-1)+1)
observed_log_cnt_n = np.log(np.sum(outs_obs,axis=-1)+1)

print(observed_log_cnts.shape)

data = {}
data["observed"] = np.squeeze(observed_log_cnt_n)
import pickle as pkl

file = open('observed.pkl', 'wb')

pkl.dump(data, file)
file.close()


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
print(pred_logits.shape)
print(pred_logits_wo_bias.shape)

pred_cnts = np.squeeze(pred_logcts)[regions_df['chr'].isin(valid_chroms)]
iso_reg_wb = IsotonicRegression().fit(pred_cnts, observed_log_cnts)

pred_cnts = np.squeeze(pred_logcts_wo_bias)[regions_df['chr'].isin(valid_chroms)]
iso_reg_wob= IsotonicRegression().fit(pred_cnts, observed_log_cnts)

data = {}
data["original"] = np.squeeze(pred_logcts)
data["observed"] = np.squeeze(observed_log_cnts)
pred_logcts = iso_reg_wb.predict(np.squeeze(pred_logcts))
data["fitted"] = np.squeeze(pred_logcts)

import pickle as pkl

file = open('important.pkl', 'wb')

pkl.dump(data, file)
file.close()

pred_logcts_wo_bias = iso_reg_wob.predict(np.squeeze(pred_logcts_wo_bias))


bigwig_helper.write_bigwig(softmax(pred_bias_logits) * (np.expand_dims(np.exp(pred_bias_logcts)[:,0],axis=1)), 
                           regions, 
                           gs, 
                           args.out_prefix + "bias.bw", 
                           args.out_prefix + "bias.stats.txt", 
                           debug_chr=args.debug_chr, 
                           use_tqdm=args.tqdm)

bigwig_helper.write_bigwig(softmax(pred_logits_wo_bias) * (np.expand_dims(np.exp(pred_logcts_wo_bias),axis=1)), 
                           regions, 
                           gs, 
                           args.out_prefix + "_wo_bias.bw", 
                           args.out_prefix + "_wo_bias.stats.txt", 
                           debug_chr=args.debug_chr, 
                           use_tqdm=args.tqdm)

bigwig_helper.write_bigwig(softmax(pred_logits) * (np.expand_dims(np.exp(pred_logcts),axis=1)),
                           regions,
                           gs,
                           args.out_prefix + "_w_bias.bw",
                           args.out_prefix + "_w_bias.stats.txt",
                           debug_chr=args.debug_chr,
                           use_tqdm=args.tqdm)
