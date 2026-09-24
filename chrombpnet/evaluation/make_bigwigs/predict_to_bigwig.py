import argparse
import pyBigWig
import numpy as np
import pandas as pd
import pyfaidx
import chrombpnet.evaluation.make_bigwigs.bigwig_helper as bigwig_helper
from chrombpnet.training.utils import model_io
import chrombpnet.training.utils.data_utils as data_utils 
import chrombpnet.training.utils.one_hot as one_hot
import h5py
import json

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def write_predictions_h5py(output_prefix, profile, logcts, coords):
    # open h5 file for writing predictions
    output_h5_fname = "{}_predictions.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("predictions")

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_center_dset =  [int(coords[i][1]) for i in range(num_examples)]

    dt = h5py.string_dtype()

    # create the "coords" group datasets
    coords_chrom_dset = coord_group.create_dataset(
        "coords_chrom", data=np.array(coords_chrom_dset, dtype=dt),
        dtype=dt, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_center", data=coords_center_dset, dtype=int, compression="gzip")

    # create the "predictions" group datasets
    profs_dset = pred_group.create_dataset(
        "profs",
        data=profile,
        dtype=float, compression="gzip")
    logcounts_dset = pred_group.create_dataset(
        "logcounts", data=logcts,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()

def compare_with_observed(bigwig, regions_df, regions, outputlen, pred_logits, pred_logcts, output_prefix):

	import chrombpnet.training.metrics as metrics 
	
	obs_bw = pyBigWig.open(bigwig)
	obs_data = data_utils.get_cts(regions_df,obs_bw,outputlen)
	
	true_counts = obs_data
	true_counts_sum = np.log(np.sum(true_counts, axis=-1)+1)
	profile_probs_predictions = softmax(pred_logits) ##
	counts_sum_predictions = np.squeeze(pred_logcts) ##
	coordinates =  [[r[0], r[-1]] for r in regions]
	
	write_predictions_h5py(output_prefix, profile_probs_predictions, counts_sum_predictions, coordinates)
	
	# store regions, their predictions and corresponding pointwise metrics
	mnll_pw, mnll_norm, jsd_pw, jsd_norm, jsd_rnd, jsd_rnd_norm, mnll_rnd, mnll_rnd_norm =  metrics.profile_metrics(true_counts,profile_probs_predictions)

	spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum, counts_sum_predictions, output_prefix, "All regions provided")
	
	metrics_dictionary={}
	metrics_dictionary["counts_metrics"] = {}
	metrics_dictionary["profile_metrics"] = {}
	metrics_dictionary["counts_metrics"]["regions"] = {}
	metrics_dictionary["counts_metrics"]["regions"]["spearmanr"] = spearman_cor
	metrics_dictionary["counts_metrics"]["regions"]["pearsonr"] = pearson_cor
	metrics_dictionary["counts_metrics"]["regions"]["mse"] = mse
	
	metrics_dictionary["profile_metrics"]["regions"] = {}
	metrics_dictionary["profile_metrics"]["regions"]["median_jsd"] = float(np.nanmedian(jsd_pw))
	metrics_dictionary["profile_metrics"]["regions"]["median_norm_jsd"] = float(np.nanmedian(jsd_norm))
	
	metrics.plot_histogram(jsd_pw, jsd_rnd, output_prefix, "All regions provided")
	
	with open(output_prefix+'_metrics.json', 'w') as fp:
		json.dump(metrics_dictionary, fp,  indent=4)


# need full paths!
def parse_args():
    parser = argparse.ArgumentParser(description="Make model predictions on given regions and output to bigwig file.Please read all parameter argument requirements. PROVIDE ABSOLUTE PATHS!")
    parser.add_argument("-bm", "--bias-model", type=str, default=None, required=False, help="Path to bias model h5")
    parser.add_argument("-cm", "--chrombpnet-model", type=str, default=None, required=False, help="Path to chrombpnet model h5")
    parser.add_argument("-cmb", "--chrombpnet-model-nb", type=str, default=None, required=False, help="Path to chrombpnet model h5")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column BED file of length = N which matches f['projected_shap']['seq'].shape[0]. The ith region in the BED file corresponds to ith entry in importance matrix. If start=2nd col, summit=10th col, then the input regions are assumed to be for [start+summit-(inputlen/2):start+summit+(inputlen/2)]. Should not be piped since it is read twice!")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
    parser.add_argument("-op", "--output-prefix", type=str, required=True, help="Output bigwig file")
    parser.add_argument("-os", "--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
    parser.add_argument("-b", "--batch-size", type=int, default=64)
    parser.add_argument("-t", "--tqdm", type=int,default=0, help="Use tqdm. If yes then you need to have it installed.")
    parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    parser.add_argument("-bw", "--bigwig", type=str, default=None, help="If provided .h5 with predictions are output along with calculated metrics considering bigwig as groundtruth.")
    args = parser.parse_args()
    assert (args.bias_model is None) + (args.chrombpnet_model is None) + (args.chrombpnet_model_nb is None) < 3, "No input model provided!"
    print(args)
    return args


def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def load_regions(args, inputlen, outputlen):
    """Regions (only those on args.debug_chr if given), their one-hot sequences, the mask of regions whose
    sequence has the full input length, and the output windows of those regions: seqs[k], regions[k] and
    regions_df[regions_used].iloc[k] are the same region."""
    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
    if args.debug_chr is not None:
        # subset before extracting sequences so that the sequences and the regions stay aligned
        debug_chr = [args.debug_chr] if isinstance(args.debug_chr, str) else list(args.debug_chr)
        regions_df = regions_df[regions_df['chr'].isin(debug_chr)].reset_index(drop=True)
        if len(regions_df) == 0:
            raise ValueError("no regions of {} are on the --debug-chr chromosomes {}".format(args.regions, debug_chr))
    print(regions_df.head())
    with pyfaidx.Fasta(args.genome) as g:
        seqs, regions_used = bigwig_helper.get_seq(regions_df, g, inputlen)
    regions = bigwig_helper.get_regions(regions_df, outputlen, regions_used) # output regions
    return regions_df, seqs, regions_used, regions

def main(args):


    if args.chrombpnet_model_nb:
        model_chrombpnet_nb = model_io.load_model_wrapper(model_h5=args.chrombpnet_model_nb)
        inputlen = int(model_chrombpnet_nb.input_shape[1])
        outputlen = int(model_chrombpnet_nb.output_shape[0][1])

        # load data
        regions_df, seqs, regions_used, regions = load_regions(args, inputlen, outputlen)
        gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
        regions_df[regions_used].to_csv(args.output_prefix + "_chrombpnet_nobias_preds.bed", sep="\t", header=False, index=False)


        pred_logits_wo_bias, pred_logcts_wo_bias = model_chrombpnet_nb.predict(seqs,
                                          batch_size = args.batch_size,
                                          verbose=True)

        pred_logits_wo_bias = np.squeeze(pred_logits_wo_bias)


        bigwig_helper.write_bigwig(softmax(pred_logits_wo_bias) * (np.expand_dims(np.exp(pred_logcts_wo_bias)[:,0],axis=1)), 
                               regions, 
                               gs, 
                               args.output_prefix + "_chrombpnet_nobias.bw", 
                               outstats_file=args.output_prefix_stats, 
                               debug_chr=args.debug_chr, 
                               use_tqdm=args.tqdm)

        if args.bigwig:
        	compare_with_observed(args.bigwig, regions_df[regions_used], regions, outputlen, 
        				pred_logits_wo_bias, pred_logcts_wo_bias, args.output_prefix+"_chrombpnet_nobias")
        	

    if args.chrombpnet_model:
        model_chrombpnet = model_io.load_model_wrapper(model_h5=args.chrombpnet_model)
        inputlen = int(model_chrombpnet.input_shape[1])
        outputlen = int(model_chrombpnet.output_shape[0][1])

        # load data
        regions_df, seqs, regions_used, regions = load_regions(args, inputlen, outputlen)
        gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)
        regions_df[regions_used].to_csv(args.output_prefix + "_chrombpnet_preds.bed", sep="\t", header=False, index=False)

        pred_logits, pred_logcts = model_chrombpnet.predict(seqs,
                                          batch_size = args.batch_size,
                                          verbose=True)


        pred_logits = np.squeeze(pred_logits)


        bigwig_helper.write_bigwig(softmax(pred_logits) * (np.expand_dims(np.exp(pred_logcts)[:,0],axis=1)),
                               regions,
                               gs,
                               args.output_prefix + "_chrombpnet.bw",
                               outstats_file=args.output_prefix_stats, 
                               debug_chr=args.debug_chr,
                               use_tqdm=args.tqdm)

        if args.bigwig:
        	compare_with_observed(args.bigwig, regions_df[regions_used], regions, outputlen, 
        				pred_logits, pred_logcts, args.output_prefix+"_chrombpnet")
        	

    if args.bias_model:
        model_bias = model_io.load_model_wrapper(model_h5=args.bias_model)
        inputlen = int(model_bias.input_shape[1])
        outputlen = int(model_bias.output_shape[0][1])

        # load data
        regions_df, seqs, regions_used, regions = load_regions(args, inputlen, outputlen)
        regions_df[regions_used].to_csv(args.output_prefix + "_bias_preds.bed", sep="\t", header=False, index=False)
        gs = bigwig_helper.read_chrom_sizes(args.chrom_sizes)

        pred_bias_logits, pred_bias_logcts = model_bias.predict(seqs,
                                          batch_size = args.batch_size,
                                          verbose=True)

        pred_bias_logits = np.squeeze(pred_bias_logits)

        bigwig_helper.write_bigwig(softmax(pred_bias_logits) * (np.expand_dims(np.exp(pred_bias_logcts)[:,0],axis=1)), 
                               regions, 
                               gs, 
                               args.output_prefix + "_bias.bw", 
                               outstats_file=args.output_prefix_stats, 
                               debug_chr=args.debug_chr, 
                               use_tqdm=args.tqdm)

        if args.bigwig:
        	compare_with_observed(args.bigwig, regions_df[regions_used], regions, outputlen, 
        				pred_bias_logits, pred_bias_logcts, args.output_prefix+"_bias")
        
    

if __name__=="__main__":
    args = parse_args()
    main(args)
