import os
import subprocess
import numpy as np
import pandas as pd
import tempfile
import util

def export_motifs(pfms, out_dir):
    """
    Exports motifs to an output directory as PFMs for MOODS.
    Arguments:
        `pfms`: a dictionary mapping keys to N x 4 NumPy arrays (N may be
            different for each PFM); `{key}.pfm` will be the name of each saved
            motif
        `out_dir`: directory to save each motif
    """
    for key, pfm in pfms.items():
        outfile = os.path.join(out_dir, "%s.pfm" % key)
        with open(outfile, "w") as f:
            for i in range(4):
                f.write(" ".join([str(x) for x in pfm[:, i]]) + "\n")


def run_moods(out_dir, reference_fasta, pval_thresh=0.01):
    """
    Runs MOODS on every `.pfm` file in `out_dir`. Outputs the results for each
    PFM into `out_dir/moods_out.csv`.
    Arguments:
        `out_dir`: directory with PFMs
        `reference_fasta`: path to reference Fasta to use
        `pval_thresh`: threshold p-value for MOODS to use
    """
    procs = []
    pfm_files = [p for p in os.listdir(out_dir) if p.endswith(".pfm")]
    comm = ["moods-dna.py"]
    comm += ["-m"]
    comm += [os.path.join(out_dir, pfm_file) for pfm_file in pfm_files]
    comm += ["-s", reference_fasta]
    comm += ["-p", str(pval_thresh)]
    comm += ["-o", os.path.join(out_dir, "moods_out.csv")]
    proc = subprocess.Popen(comm)
    proc.wait()


def moods_hits_to_bed(moods_out_csv_path, moods_out_bed_path):
    """
    Converts MOODS hits into BED file.
    """
    f = open(moods_out_csv_path, "r")
    g = open(moods_out_bed_path, "w")
    warn = True
    for line in f:
        tokens = line.split(",")
        try:
            # The length of the interval is the length of the motif
            g.write("\t".join([
                tokens[0].split()[0], tokens[2],
                str(int(tokens[2]) + len(tokens[5])), tokens[1][:-4], tokens[3],
                tokens[4]
            ]) + "\n")
        except ValueError:
            # If a line is formatted incorrectly, skip it and warn once
            if warn:
                print("Found bad line: " + line)
                warn = False
            pass
        # Note: depending on the Fasta file and version of MOODS, only keep the
        # first token of the "chromosome"
    f.close()
    g.close()


def filter_hits_for_peaks(
    moods_out_bed_path, filtered_hits_path, peak_bed_path
):
    """
    Filters MOODS hits for only those that (fully) overlap a particular set of
    peaks. `peak_bed_path` must be a BED file; only the first 3 columns are
    used. A new column is added to the resulting hits: the index of the peak in
    `peak_bed_path`. If `peak_bed_path` has repeats, the later index is kept.
    """
    # First filter using bedtools intersect, keeping track of matches
    temp_file = filtered_hits_path + ".tmp"
    comm = ["bedtools", "intersect"]
    comm += ["-wa", "-wb"]
    comm += ["-f", "1"]  # Require the entire hit to overlap with peak
    comm += ["-a", moods_out_bed_path]
    comm += ["-b", peak_bed_path]
    with open(temp_file, "w") as f:
        proc = subprocess.Popen(comm, stdout=f)
        proc.wait()

    # Create mapping of peaks to indices in `peak_bed_path`
    peak_table = pd.read_csv(
        peak_bed_path, sep="\t", header=None, index_col=False,
        usecols=[0, 1, 2], names=["chrom", "start", "end"]
    )
    peak_keys = (
        peak_table["chrom"] + ":" + peak_table["start"].astype(str) + "-" + \
        peak_table["end"].astype(str)
    ).values
    peak_index_map = {k : str(i) for i, k in enumerate(peak_keys)}

    # Convert last three columns to peak index
    f = open(temp_file, "r")
    g = open(filtered_hits_path, "w")
    for line in f:
        tokens = line.strip().split("\t")
        g.write("\t".join((tokens[:-3])))
        peak_index = peak_index_map["%s:%s-%s" % tuple(tokens[-3:])]
        g.write("\t" + peak_index + "\n")
    f.close()
    g.close()


def compute_hits_importance_scores(
    hits_bed_path, shap_scores_hdf5_path, peak_table, input_len, peak_bed_path, out_path, method
):
    """
    For each MOODS hit, computes the hit's importance score as the ratio of the
    hit's average importance score to the total importance of the sequence.
    Arguments:
        `hits_bed_path`: path to BED file output by `collapse_hits`
            without the p-value column
        `shap_scores_hdf5_path`: an HDF5 of DeepSHAP scores of peak regions
            measuring importance
        `peak_bed_path`: BED file of peaks; we require that these coordinates
            must match the DeepSHAP score coordinates exactly
        `out_path`: path to output the resulting table
    Each of the DeepSHAP score HDF5s must be of the form:
        `coords_chrom`: N-array of chromosome (string)
        `coords_start`: N-array
        `coords_end`: N-array
        `hyp_scores`: N x L x 4 array of hypothetical importance scores
        `input_seqs`: N x L x 4 array of one-hot encoded input sequences
    Outputs an identical hit BED with an extra column for the importance score
    fraction.
    """
    #_, imp_scores, _, coords = util.import_shap_scores_custom(
    #    shap_scores_hdf5_path, peak_table_in, center_cut_size=input_len
    #)
    
    imp_scores = util.import_shap_scores_from_bigwig(
        shap_scores_hdf5_path, peak_table, input_len
    )
    
    #print(imp_scores.shape)
    peak_table = pd.read_csv(
        peak_bed_path, sep="\t", header=None, index_col=False,
        usecols=[0, 1, 2], names=["peak_chrom", "peak_start", "peak_end"]
    )
    #assert np.all(coords == peak_table.values)
    
    #print(peak_table.head())

    hit_table = pd.read_csv(
        hits_bed_path, sep="\t", header=None, index_col=False,
        names=["chrom", "start", "end", "key", "strand", "score", "peak_index"]
    )

    #print(hit_table.head())
    
    # Merge in the peak starts/ends to the hit table
    merged_hits = pd.merge(
        hit_table, peak_table, left_on="peak_index", right_index=True
    )

    #print(merged_hits.head())

    # Important! Reset the indices of `merged_hits` after merging, otherwise
    # iteration over the rows won't be in order
    merged_hits = merged_hits.reset_index(drop=True)

    # Compute start and end of each motif relative to the peak
    merged_hits["motif_rel_start"] = \
        merged_hits["start"] - merged_hits["peak_start"]
    merged_hits["motif_rel_end"] = \
        merged_hits["end"] - merged_hits["peak_start"]

    # Careful! Because of the merging step that only kept the top peak hit, some
    # hits might overrun the edge of the peak; we limit the motif hit indices
    # here so they stay in the peak; this should not be a common occurrence
    merged_hits["peak_min"] = 0
    merged_hits["peak_max"] = \
        merged_hits["peak_end"] - merged_hits["peak_start"]
    merged_hits["motif_rel_start"] = \
        merged_hits[["motif_rel_start", "peak_min"]].max(axis=1)
    merged_hits["motif_rel_end"] = \
        merged_hits[["motif_rel_end", "peak_max"]].min(axis=1)
    del merged_hits["peak_min"]
    del merged_hits["peak_max"]

    # Get score of each motif hit as average importance over the hit, divided
    # by the total score of the sequence
    scores = np.empty(len(merged_hits))
    for peak_index, group in merged_hits.groupby("peak_index"):
        # Iterate over grouped table by peak
        score_track = np.abs(imp_scores[peak_index])
        total_score = np.mean(score_track)
        for i, row in group.iterrows():
            if method=="norm":
                scores[i] = np.sum(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                )/total_score
            if method=="sum":
                scores[i] = np.sum(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                )
            if method=="mean_norm":
                scores[i] = np.mean(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                )/total_score
            if method=="mean_sum":
                scores[i] = np.mean(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                )
    merged_hits["imp_frac_score"] = scores
    new_hit_table = merged_hits[[
        "chrom", "start", "end", "key", "strand", "score", "peak_index",
        "imp_frac_score"
    ]]
    new_hit_table.to_csv(out_path, sep="\t", header=False, index=False)


def import_moods_hits(hits_bed):
    """
    Imports the MOODS hits as a single Pandas DataFrame.
    Returns a Pandas DataFrame with the columns: chrom, start, end, key, strand,
    score, peak_index, imp_frac_score
    `key` is the name of the originating PFM, and `length` is its length.
    """
    hit_table = pd.read_csv(
        hits_bed, sep="\t", header=None, index_col=False,
        names=[
            "chrom", "start", "end", "key", "strand", "score", "peak_index",
            "imp_frac_score"
        ]
    )
    return hit_table


def get_moods_hits(
    pfm_dict, reference_fasta, peak_bed_path, shap_scores_hdf5_path, peak_table_in,
    expand_peak_length=None, moods_pval_thresh=0.01, temp_dir=None, method="sum"
):
    """
    From a dictionary of PFMs, runs MOODS and returns the result as a Pandas
    DataFrame.
    Arguments:
        `pfm_dict`: a dictionary mapping keys to N x 4 NumPy arrays (N may be
            different for each PFM); the key will be the name of each motif
        `reference_fasta`: path to reference Fasta to use
        `peak_bed_path`: path to peaks BED file; only keeps MOODS hits from
            these intervals; must be in NarrowPeak format
        `shap_scores_hdf5_path`: an HDF5 of DeepSHAP scores of peak regions
            measuring importance
        `peak_bed_path`: BED file of peaks; we require that these coordinates
            must match the DeepSHAP score coordinates exactly
        `expand_peak_length`: if given, expand the peaks (centered at summits)
            to this length
        `moods_pval_thresh`: threshold p-value for MOODS to use
        `temp_dir`: a temporary directory to store intermediates; defaults to
            a randomly created directory
    Each of the DeepSHAP score HDF5s must be of the form:
        `
        _chrom`: N-array of chromosome (string)
        `coords_start`: N-array
	`coords_end`: N-array
	`hyp_scores`: N x L x 4 array of hypothetical importance scores
	`input_seqs`: N x L x 4 array of one-hot encoded input sequences
    The coordinates of the DeepSHAP scores must be identical, and must match
    the peaks in the BED file (after expansion, if specified).
    """
    if temp_dir is None:
        temp_dir_obj = tempfile.TemporaryDirectory()
        temp_dir = temp_dir_obj.name
    else:
        os.makedirs(temp_dir, exist_ok=True)
        temp_dir_obj = None

    pfm_keys = list(pfm_dict.keys())
    
    # Create PFM files
    export_motifs(pfm_dict, temp_dir)

    print("Running moods...")
    # Run MOODS
    run_moods(temp_dir, reference_fasta, pval_thresh=moods_pval_thresh)

    # Convert MOODS output into BED file
    moods_hits_to_bed(
        os.path.join(temp_dir, "moods_out.csv"),
        os.path.join(temp_dir, "moods_out.bed")
    )

    print("Filtering hits that overlap with peaks...")
    # Filter hits for those that overlap peaks
    filter_hits_for_peaks(
        os.path.join(temp_dir, "moods_out.bed"),
        os.path.join(temp_dir, "moods_filtered.bed"),
        peak_bed_path
    )

    
    print("Compute importance scores in hits")
    compute_hits_importance_scores(
        os.path.join(temp_dir, "moods_filtered.bed"),
        shap_scores_hdf5_path, peak_table_in, expand_peak_length, peak_bed_path,
        os.path.join(temp_dir, "moods_filtered_scored.bed"), method
    )

    hit_table = import_moods_hits(
        os.path.join(temp_dir, "moods_filtered_scored.bed")
    )

    if temp_dir_obj is not None:
        temp_dir_obj.cleanup()

    return hit_table


def get_moods_hits_from_filtered(
    pfm_dict, reference_fasta, peak_bed_path, shap_scores_hdf5_path, peak_table_in,
    expand_peak_length=None, moods_pval_thresh=0.001, temp_dir=None, method="sum"
):
    """
    From a dictionary of PFMs, runs MOODS and returns the result as a Pandas
    DataFrame.
    Arguments:
        `pfm_dict`: a dictionary mapping keys to N x 4 NumPy arrays (N may be
            different for each PFM); the key will be the name of each motif
        `reference_fasta`: path to reference Fasta to use
        `peak_bed_path`: path to peaks BED file; only keeps MOODS hits from
            these intervals; must be in NarrowPeak format
        `shap_scores_hdf5_path`: an HDF5 of DeepSHAP scores of peak regions
            measuring importance
        `peak_bed_path`: BED file of peaks; we require that these coordinates
            must match the DeepSHAP score coordinates exactly
        `expand_peak_length`: if given, expand the peaks (centered at summits)
            to this length
        `moods_pval_thresh`: threshold p-value for MOODS to use
        `temp_dir`: a temporary directory to store intermediates; defaults to
            a randomly created directory
    Each of the DeepSHAP score HDF5s must be of the form:
        `coords_chrom`: N-array of chromosome (string)
        `coords_start`: N-array
	`coords_end`: N-array
	`hyp_scores`: N x L x 4 array of hypothetical importance scores
	`input_seqs`: N x L x 4 array of one-hot encoded input sequences
    The coordinates of the DeepSHAP scores must be identical, and must match
    the peaks in the BED file (after expansion, if specified).
    """


    print("calculate scores for hits...")
    compute_hits_importance_scores(
        os.path.join(temp_dir, "moods_filtered.bed"),
        shap_scores_hdf5_path, peak_table_in, expand_peak_length, peak_bed_path,
        os.path.join(temp_dir, "moods_filtered_scored.bed"), method
    )

    hit_table = import_moods_hits(
        os.path.join(temp_dir, "moods_filtered_scored.bed")
    )

 
    return hit_table
