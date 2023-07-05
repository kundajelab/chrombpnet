import deepdish as dd
import json
import numpy as np
import tensorflow as tf
import pandas as pd
import shap
from tensorflow import keras
import pyfaidx
import shutil
import errno
import os
from utils import argmanager, data_utils
from utils.loss import multinomial_nll
from utils.shap_utils import *

# disable eager execution so shap deep explainer wont break
tf.compat.v1.disable_eager_execution()


def generate_shap_dict(seqs, scores):
    assert(seqs.shape==scores.shape)
    assert(seqs.shape[2]==4)

    # construct a dictionary for the raw shap scores and the
    # the projected shap scores
    # MODISCO workflow expects one hot sequences with shape (None,4,inputlen)
    d = {
            'raw': {'seq': np.transpose(seqs, (0, 2, 1))},
            'shap': {'seq': np.transpose(scores, (0, 2, 1))},
            'projected_shap': {'seq': np.transpose(seqs*scores, (0, 2, 1))}
        }

    return d

def interpret(model, seqs, output_prefix, profile_or_counts):
    print("Seqs dimension : {}".format(seqs.shape))

    is_chrombpnet_arch = isinstance(model.input, list)
    outlen = model.output_shape[0][1]

    if is_chrombpnet_arch:
        assert(len(model.input)==3)
        # 3 inputs [seq, bias_logits (None, seqlen), bias_cts (None,1)]
        profile_model_input = [model.input[0], model.input[1]]
        profile_input = [seqs, np.zeros((seqs.shape[0], outlen))]
        counts_model_input = [model.input[0], model.input[2]]
        counts_input = [seqs, np.zeros((seqs.shape[0], 1))]
    else:
        # is bias arch (only seq input)
        profile_model_input = model.input
        profile_input = seqs
        counts_model_input = model.input
        counts_input = seqs

    if "counts" in profile_or_counts:
        profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
            (counts_model_input, tf.reduce_sum(model.outputs[1], axis=-1)),
            shuffle_several_times,
            combine_mult_and_diffref=combine_mult_and_diffref)

        print("Generating 'counts' shap scores")
        counts_shap_scores = profile_model_counts_explainer.shap_values(
            counts_input, progress_message=100)

        if is_chrombpnet_arch:
            counts_shap_scores = counts_shap_scores[0]

        counts_scores_dict = generate_shap_dict(seqs, counts_shap_scores)

        # save the dictionary in HDF5 formnat
        print("Saving 'counts' scores")
        dd.io.save("{}.counts_scores.h5".format(output_prefix),
                    counts_scores_dict,
                    compression='blosc')

        del counts_shap_scores, counts_scores_dict

    if "profile" in profile_or_counts:
        weightedsum_meannormed_logits = get_weightedsum_meannormed_logits(model)
        profile_model_profile_explainer = shap.explainers.deep.TFDeepExplainer(
            (profile_model_input, weightedsum_meannormed_logits),
            shuffle_several_times,
            combine_mult_and_diffref=combine_mult_and_diffref)

        print("Generating 'profile' shap scores")
        profile_shap_scores = profile_model_profile_explainer.shap_values(
            profile_input, progress_message=100)

        if is_chrombpnet_arch:
            # take scores corresponding to sequence
            profile_shap_scores = profile_shap_scores[0]

        profile_scores_dict = generate_shap_dict(seqs, profile_shap_scores)

        # save the dictionary in HDF5 formnat
        print("Saving 'profile' scores")
        dd.io.save("{}.profile_scores.h5".format(output_prefix),
                    profile_scores_dict,
                    compression='blosc')


def main():
    # parse the command line arguments
    args = argmanager.fetch_interpret_args()

    # check if the output directory exists
    if not os.path.exists(os.path.dirname(args.output_prefix)):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.dirname(args.output_prefix))

    # write all the command line arguments to a json file
    with open("{}.args.json".format(args.output_prefix), "w") as fp:
        json.dump(vars(args), fp, ensure_ascii=False, indent=4)

    regions_df = pd.read_csv(args.regions, sep='\t', names=data_utils.NARROWPEAK_SCHEMA)

    if args.debug_chr:
        regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]
        regions_df.to_csv("{}.interpreted_regions.bed".format(args.output_prefix), header=False, sep='\t')
    else:
        # copy regions bed to output directory
        shutil.copy(args.regions, "{}.interpreted_regions.bed".format(args.output_prefix))

    # load the model
    with keras.utils.CustomObjectScope({'multinomial_nll':multinomial_nll, 'tf':tf}):
        model = keras.models.load_model(args.model)

    # infer input length
    if isinstance(model.input, list):
        inputlen = model.input_shape[0][1] # if chrombpnet model (3 inputs, first is sequence)
    else:
        inputlen = model.input_shape[1] # if bias model (1 input only)

    # load sequences
    # NOTE: it will pull out sequences of length inputlen
    #       centered at the summit (start + 10th column)
    genome = pyfaidx.Fasta(args.genome)
    seqs = data_utils.get_seq(regions_df, genome, inputlen)
    genome.close()

    interpret(model, seqs, args.output_prefix, args.profile_or_counts)

if __name__ == '__main__':
    main()
