import h5py
import numpy as np
import pandas as pd
import modisco
import tqdm
import io
import base64
import urllib
import vdom.helpers as vdomh
import deepdish
import pyBigWig

def import_shap_scores_from_bigwig(
    shap_scores_bw, peak_table, input_len
):
    bw=pyBigWig.open(shap_scores_bw)
    scores=[]
    for i,r in peak_table.iterrows():
        val = np.nan_to_num(bw.values(r["chrom"],r["peak_start"],r["peak_end"]))
        scores.append(val)
    return np.array(scores) 
       

def import_profiles(preds_path):
    """
    Imports the set of profile predictions.
    Arguments:
        `preds_path`: path to predictions/performance metrics of the model
    Returns an M x T x O x 2 array of true profile counts, an M x T x O x 2
    array of predicted profile probabilities, and an M x 3 object array of
    corresponding coordinates.
    """
    with h5py.File(preds_path, "r") as f:
        num_seqs, num_tasks, input_length, _ = f["predictions"]["true_profs"].shape
        batch_size = min(1000, num_seqs)
        num_batches = int(np.ceil(num_seqs / batch_size))
        
        true_profs = np.empty((num_seqs, num_tasks, input_length, 2))
        pred_profs = np.empty((num_seqs, num_tasks, input_length, 2))
        coords = np.empty((num_seqs, 3), dtype=object)
        
        for i in tqdm.notebook.trange(num_batches, desc="Importing predictions"):
            batch_slice = slice(i * batch_size, (i + 1) * batch_size)
            true_profs[batch_slice] = f["predictions"]["true_profs"][batch_slice]
            pred_profs[batch_slice] = np.exp(f["predictions"]["log_pred_profs"][batch_slice])
            coords[batch_slice, 0] = f["coords"]["coords_chrom"][batch_slice].astype(str)
            coords[batch_slice, 1] = f["coords"]["coords_start"][batch_slice]
            coords[batch_slice, 2] = f["coords"]["coords_end"][batch_slice]
    
    return true_profs, pred_profs, coords


def import_tfmodisco_results(tfm_results_path, hyp_scores, one_hot_seqs, input_center_cut_size):
    """
    Imports the TF-MoDISco results object.
    Arguments:
        `tfm_results_path`: path to HDF5 containing TF-MoDISco results
        `hyp_scores`: hypothetical importance scores used for this run
        `one_hot_seqs`: input sequences used for this run
        `input_center_cut_size`: centered cut size of SHAP scores used
    """ 
    # Everything should already be cut to `input_center_cut_size`
    act_scores = hyp_scores * one_hot_seqs
    
    track_set = modisco.tfmodisco_workflow.workflow.prep_track_set(
        task_names=["task0"],
        contrib_scores={"task0": act_scores},
        hypothetical_contribs={"task0": hyp_scores},
        one_hot=one_hot_seqs
    )
    
    with h5py.File(tfm_results_path,"r") as f:
        return modisco.tfmodisco_workflow.workflow.TfModiscoResults.from_hdf5(f, track_set=track_set)


def import_peak_table_custom(peak_bed_paths):
    tables = []
    for peak_bed_path in peak_bed_paths:
        table = pd.read_csv(
            peak_bed_path, sep="\t", header=None,  # Infer compression
            names=[
                "chrom", "peak_start", "peak_end"]
        )
        tables.append(table)
    return pd.concat(tables)

def import_peak_table(peak_bed_paths):
    tables = []
    for peak_bed_path in peak_bed_paths:
        table = pd.read_csv(
            peak_bed_path, sep="\t", header=None,  # Infer compression
            names=[
                "chrom", "peak_start", "peak_end", "name", "score",
                "strand", "signal", "pval", "qval", "summit_offset"]
        )
        # Add summit location column
        table["summit"] = table["peak_start"] + table["summit_offset"]
        tables.append(table)
    return pd.concat(tables)


BACKGROUND_FREQS = np.array([0.25, 0.25, 0.25, 0.25])
def info_content(track, pseudocount=0.001):
    """
    Given an L x 4 track, computes information content for each base and
    returns it as an L-array.
    """
    num_bases = track.shape[1]
    # Normalize track to probabilities along base axis
    track_norm = (track + pseudocount) / (np.sum(track, axis=1, keepdims=True) + (num_bases * pseudocount))
    ic = track_norm * np.log2(track_norm / np.expand_dims(BACKGROUND_FREQS, axis=0))
    return np.sum(ic, axis=1)


def pfm_to_pwm(pfm):
    ic = info_content(pfm)
    return pfm * np.expand_dims(ic, axis=1)


def trim_motif(pfm, motif, min_ic=0.2, pad=0):
    """
    Given the PFM and motif (both L x 4 arrays) (the motif could be the
    PFM itself), trims `motif` by cutting off flanks of low information
    content in `pfm`. `min_ic` is the minimum required information
    content. If specified this trimmed motif will be extended on either
    side by `pad` bases.
    If no base passes the `min_ic` threshold, then no trimming is done.
    """
    # Trim motif based on information content
    ic = info_content(pfm)
    pass_inds = np.where(ic >= min_ic)[0]  # Cut off flanks with less than min_ic IC
    
    if not pass_inds.size:
        return motif

    # Expand trimming to +/- pad bp on either side
    start, end = max(0, np.min(pass_inds) - pad), min(len(pfm), np.max(pass_inds) + pad + 1)
    return motif[start:end]

def trim_motif_new(cwm, motif, trim_threshold=0.3):
    """
    Given the PFM and motif (both L x 4 arrays) (the motif could be the
    PFM itself), trims `motif` by cutting off flanks of low information
    content in `pfm`. `min_ic` is the minimum required information
    content. If specified this trimmed motif will be extended on either
    side by `pad` bases.
    If no base passes the `min_ic` threshold, then no trimming is done.
    """
    
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
    pass_inds = np.where(score >= trim_thresh)[0]
    trimmed = motif[np.min(pass_inds): np.max(pass_inds) + 1]
 
    if not trimmed.size:
        return motif
    
    return trimmed


def figure_to_vdom_image(figure):
    buf = io.BytesIO()
    figure.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
        
    return vdomh.div(
        vdomh.img(src='data:image/png;base64,' + urllib.parse.quote(string)),
        style={"display": "inline-block"}
    )
