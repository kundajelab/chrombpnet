# Adapted from chrombpnet-lite

import deepdish as dd
import json
import numpy as np
import tensorflow as tf
import pandas as pd
import shap
import pyfaidx
import shutil
import errno
import os
import argparse
import chrombpnet.evaluation.interpret.shap_utils as shap_utils
import chrombpnet.evaluation.interpret.input_utils as input_utils

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

# disable eager execution so shap deep explainer wont break
tf.compat.v1.disable_eager_execution()

def fetch_interpret_args():
    parser = argparse.ArgumentParser(description="get sequence contribution scores for the model")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-o", "--output-prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-d", "--debug_chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    parser.add_argument("-p", "--profile_or_counts", nargs="+", type=str, default=["counts", "profile"], choices=["counts", "profile"],
                        help="use either counts or profile or both for running shap")

    args = parser.parse_args()
    return args


def generate_shap_dict(seqs, scores):
    assert(seqs.shape==scores.shape)
    assert(seqs.shape[2]==4)

    # construct a dictionary for the raw shap scores and the
    # the projected shap scores
    # MODISCO workflow expects one hot sequences with shape (None,4,inputlen)
    d = {
            'raw': {'seq': np.transpose(seqs, (0, 2, 1)).astype(np.int8)},
            'shap': {'seq': np.transpose(scores, (0, 2, 1)).astype(np.float16)},
            'projected_shap': {'seq': np.transpose(seqs*scores, (0, 2, 1)).astype(np.float16)}
        }

    return d

def interpret(model, seqs, output_prefix, profile_or_counts):
    print("Seqs dimension : {}".format(seqs.shape))

    outlen = model.output_shape[0][1]

    profile_model_input = model.input
    profile_input = seqs
    counts_model_input = model.input
    counts_input = seqs

    if "counts" in profile_or_counts:
        profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
            (counts_model_input, tf.reduce_sum(model.outputs[1], axis=-1)),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)

        print("Generating 'counts' shap scores")
        counts_shap_scores = profile_model_counts_explainer.shap_values(
            counts_input, progress_message=100)

        counts_scores_dict = generate_shap_dict(seqs, counts_shap_scores)

        # save the dictionary in HDF5 formnat
        print("Saving 'counts' scores")
        dd.io.save("{}.counts_scores.h5".format(output_prefix),
                    counts_scores_dict,
                    compression='blosc')

        del counts_shap_scores, counts_scores_dict

    if "profile" in profile_or_counts:
        weightedsum_meannormed_logits = shap_utils.get_weightedsum_meannormed_logits(model)
        profile_model_profile_explainer = shap.explainers.deep.TFDeepExplainer(
            (profile_model_input, weightedsum_meannormed_logits),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)

        print("Generating 'profile' shap scores")
        profile_shap_scores = profile_model_profile_explainer.shap_values(
            profile_input, progress_message=100)

        profile_scores_dict = generate_shap_dict(seqs, profile_shap_scores)

        # save the dictionary in HDF5 formnat
        print("Saving 'profile' scores")
        dd.io.save("{}.profile_scores.h5".format(output_prefix),
                    profile_scores_dict,
                    compression='blosc')


def main(args):

    # check if the output directory exists
    #if not os.path.exists(os.path.dirname(args.output_prefix)):
    #    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.dirname(args.output_prefix))

    # write all the command line arguments to a json file
    with open("{}.interpret.args.json".format(args.output_prefix), "w") as fp:
        json.dump(vars(args), fp, ensure_ascii=False, indent=4)

    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)

    if args.debug_chr:
        regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
    
    model = input_utils.load_model_wrapper(args)
 
    # infer input length
    inputlen = model.input_shape[1] # if bias model (1 input only)
    print("inferred model inputlen: ", inputlen)

    # load sequences
    # NOTE: it will pull out sequences of length inputlen
    #       centered at the summit (start + 10th column) and peaks used after filtering

    genome = pyfaidx.Fasta(args.genome)
    seqs, peaks_used = input_utils.get_seq(regions_df, genome, inputlen)
    genome.close()

    regions_df[peaks_used].to_csv("{}.interpreted_regions.bed".format(args.output_prefix), header=False, index=False, sep='\t')

    interpret(model, seqs, args.output_prefix, args.profile_or_counts)

if __name__ == '__main__':
    # parse the command line arguments
    args = fetch_interpret_args()
    main(args)

