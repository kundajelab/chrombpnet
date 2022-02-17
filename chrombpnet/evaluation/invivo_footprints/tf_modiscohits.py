import os
import numpy as np
import pandas as pd
from modisco.hit_scoring import densityadapted_hitscoring
from modisco.util import compute_per_position_ic
import chrombpnet.evaluation.invivo_footprints.run_tfmodisco as run_tfmodisco
import click

def import_tfmodisco_hits(hits_bed):
    """
    Imports the TF-MoDISco hits as a single Pandas DataFrame.
    The `key` column is the name of the originating PFM, and `peak_index` is the
    index of the peak file from which it was originally found.
    """
    hit_table = pd.read_csv(
        hits_bed, sep="\t", header=None, index_col=False,
        names=[
            "chrom", "start", "end", "key", "strand", "peak_index",
            "imp_total_signed_score", "imp_total_score", "imp_frac_score",
            "imp_ic_avg_score", "agg_sim", "mod_delta", "mod_precision",
            "mod_percentile", "fann_perclasssum_perc", "fann_perclassavg_perc"
        ]
    )
    return hit_table


@click.command()
@click.option(
    "-o", "--outdir", required=True, help="Path to output directory"
)

@click.option(
    "-i", "--input-length", default=2114,
    help="Length of input sequences for importance scores"
)
@click.option(
    "-c", "--center-cut-size", default=1000,
    help="Length of sequence that was used to run TF-MoDISco"
)
@click.option(
    "--keep-non-acgt", is_flag=True,
    help="If given, don't remove non-ACGT score tracks"
)
@click.option(
    "-m", "--min-ic", default=0.2,
    help="Information content cut-off to use to trim motif hits"
)
@click.option(
    "-mc", "--metacluster-ind", default=0,
    help="Index of the metacluster whose patterns to use for motif assignment; defaults to metacluster 0"
)
@click.option(
    "-p", "--pattern-inds", default=None, type=str,
    help="Comma-delimited list of pattern indices in the metacluster to use for motif assignment; defaults to all patterns in the metacluster"
)
@click.argument("shap_scores_path", nargs=1)
@click.argument("tfm_results_path", nargs=1)
@click.argument("peak_bed_path", nargs=1)
def main(
    shap_scores_path, tfm_results_path, peak_bed_path, outdir,
    input_length, center_cut_size, keep_non_acgt, min_ic, metacluster_ind,
    pattern_inds
):
    assert metacluster_ind in (0, 1)
    if pattern_inds is not None:
        pattern_inds = [int(x) for x in pattern_inds.split(",")]

    os.makedirs(outdir, exist_ok=True)

    # Import peaks
    peak_table = pd.read_csv(
        peak_bed_path, sep="\t", header=None, index_col=False,
        usecols=[0, 1, 2, 9],
        names=["peak_chrom", "peak_start", "peak_end", "summit_offset"]
    )
    
    # Expand peaks to input length
    peak_table["peak_start"] = \
        (peak_table["peak_start"] + peak_table["summit_offset"]) - \
        (input_length // 2)
    peak_table["peak_end"] = peak_table["peak_start"] + input_length
    
    print("Importing DeepSHAP scores and TF-MoDISco results...")
    hyp_scores, act_scores, one_hot_seqs, imp_coords = \
        run_tfmodisco.import_shap_scores_part2(
            shap_scores_path, peak_table, center_cut_size=center_cut_size
    )
    tfm_results = run_tfmodisco.import_tfmodisco_results(
        tfm_results_path, hyp_scores, one_hot_seqs, center_cut_size
    )

    assert np.all(imp_coords[:, 2] - imp_coords[:, 1] == input_length)

    peak_table = peak_table.reset_index().drop_duplicates(
        ["peak_chrom", "peak_start", "peak_end"]
    )
    # Importantly, we add the index column before dropping duplicates
   
    print("Matching up DeepSHAP coordinates and peak coordinates...")
    imp_coords_table = pd.DataFrame(
        imp_coords, columns=["chrom", "start", "end"]
    ).reset_index().drop_duplicates(["chrom", "start", "end"])
    # Importantly, we add the index column before dropping duplicates

    # Map peak indices to importance score tracks
    matched_inds = peak_table.merge(
        imp_coords_table, how="inner", 
        # Inner join: can't call hits if there's no importance score track,
        # and don't bother if it's not a peak
        left_on=["peak_chrom", "peak_start", "peak_end"],
        right_on=["chrom", "start", "end"]
    )[["index_x", "index_y"]].values
    
    # `matched_inds` is an N x 2 array, where each pair is
    # (peak index, score index)
    # Sort by score index
    matched_inds = matched_inds[np.argsort(matched_inds[:, 1])]

    # Limit the importance scores to only those which matched to a peak
    score_inds = matched_inds[:, 1]
    hyp_scores_matched = hyp_scores[score_inds]
    act_scores_matched = act_scores[score_inds]
    one_hot_seqs_matched = one_hot_seqs[score_inds]
    
    example_to_peak_index = matched_inds[:, 0]
    # `example_to_peak_index` is an array such that if `i` is the index of
    # a sequence in `*_scores_matched`, then `example_to_peak_index[i]` is the
    # index of the matching peak
   
    print("Preparing the hit scorer...")
    # Only do the first metacluster (positive scores)
    patterns = tfm_results.metacluster_idx_to_submetacluster_results[
        "metacluster_%d" % metacluster_ind
    ].seqlets_to_patterns_result.patterns

    # If specified, use only specific patterns in the metacluster
    if pattern_inds is None:
        pattern_inds = list(range(len(patterns)))
    else:
        patterns = [patterns[i] for i in pattern_inds]

    bg_freq = np.mean(one_hot_seqs_matched, axis=(0, 1))
        
    # Verify that every pattern has sufficiently high IC
    for pattern in patterns:
        pfm = pattern["sequence"].fwd
        ic = compute_per_position_ic(pfm, bg_freq, 0.001)
        assert np.sum(ic >= min_ic) > 0, "The given IC threshold results in an empty motif"

    # Instantiate the hit scorer
    hit_scorer = densityadapted_hitscoring.MakeHitScorer(
        patterns=patterns,
        target_seqlet_size=25,
        bg_freq=bg_freq,
        task_names_and_signs=[("task0", 1 if metacluster_ind == 0 else -1)],
        n_cores=10,
        additional_seqletscorer_kwargs={"ic_trim_threshold": min_ic}
    )
    
    # Set seqlet identification method
    hit_scorer.set_coordproducer(
        contrib_scores={"task0": act_scores_matched},
        core_sliding_window_size=5,
        target_fdr=0.2,
        min_passing_windows_frac=0.03,
        max_passing_windows_frac=0.2,
        separate_pos_neg_thresholds=False,                             
        max_seqlets_total=np.inf
    )

    # Map pattern index to motif key
    motif_keys = ["%d_%d" % (metacluster_ind, i) for i in pattern_inds]

    print("Starting hit scoring...")
    batch_size = 1024
    num_batches = int(np.ceil(len(act_scores_matched) / batch_size))
    rows = []

    for i in range(num_batches):
        print("\tScoring batch %d/%d" % (i + 1, num_batches))
        batch_slice = slice(i * batch_size, (i + 1) * batch_size)
        example_to_matches, pattern_to_matches = hit_scorer(
            contrib_scores={"task0": act_scores_matched[batch_slice]},
            hypothetical_contribs={"task0": hyp_scores_matched[batch_slice]},
            one_hot=one_hot_seqs_matched[batch_slice],
            hits_to_return_per_seqlet=1
        )
        
        offset = i * batch_size
        for example_index, match_list in example_to_matches.items():
            for match in match_list:
                rows.append([
                    match.exampleidx + offset, match.patternidx, match.start,
                    match.end, match.is_revcomp, match.aggregate_sim,
                    match.mod_delta, match.mod_precision, match.mod_percentile,
                    match.fann_perclasssum_perc, match.fann_perclassavg_perc
                ])
    
    # Collate the matches together into a big table
    colnames = [
        "example_index", "pattern_index", "start", "end", "revcomp",
        "agg_sim", "mod_delta", "mod_precision", "mod_percentile",
        "fann_perclasssum_perc", "fann_perclassavg_perc"
    ]
    match_table = pd.DataFrame(rows, columns=colnames)

    # Save raw table
    match_table.to_csv(
        os.path.join(outdir, "tfm_matches_raw.tsv"), sep="\t", header=True,
        index=False
    )

    print("Cleaning up matches...")
    # Convert example index to peak index
    match_table["peak_index"] = example_to_peak_index[
        match_table["example_index"]
    ]
    
    # Convert pattern index to motif key
    match_table["key"] = np.array(motif_keys)[match_table["pattern_index"]]
    
    # Convert revcomp to strand
    # Note we are assuming that the input scores were all positive strand
    match_table["strand"] = match_table["revcomp"].map({True: "-", False: "+"})

    # Save the start/end as other columns, which match the score coordinates
    match_table["score_start"] = match_table["start"]
    match_table["score_end"] = match_table["end"]
    
    # Convert start/end of motif hit to genomic coordinate
    # `peak_starts[i] == j` is such that if `i` is a peak index, `j` is the peak
    # start in genomic coordinate space
    peak_starts = np.empty(np.max(peak_table["index"]) + 1, dtype=int)
    peak_starts[peak_table["index"]] = peak_table["peak_start"]
    # Now reduce `peak_starts` to match `match_table` exactly
    peak_starts = peak_starts[match_table["peak_index"]]
    offset = (input_length - center_cut_size) // 2

    match_table["chrom"] = peak_table["peak_chrom"].loc[
        match_table["peak_index"]
    ].reset_index(drop=True)
    # Note: "peak_chrom" was an index column so we need to drop that before
    # setting it as a value
    match_table["start"] = match_table["start"] + offset + peak_starts
    match_table["end"] = match_table["end"] + offset + peak_starts

    # Trim each motif hit to be only the size of the core motif (determined by
    # IC)
    ic_dict = {}  # Save (trimmed) IC for each motif
    hit_patterns = hit_scorer.trimmed_subclustered_patterns
    for i, pattern in enumerate(hit_patterns):
        motif_key = motif_keys[i]
        pfm = pattern["sequence"].fwd
        ic = compute_per_position_ic(pfm, bg_freq, 0.001)
        pass_inds = np.where(ic >= min_ic)[0]
        if not pass_inds.size:
            continue

        start, end = np.min(pass_inds), np.max(pass_inds) + 1
        length = end - start
        rc_start = len(pfm) - end

        ic_dict[motif_key] = ic[start:end]

        motif_mask = match_table["key"] == motif_key
        pos_mask = match_table["strand"] == "+"
        match_table.loc[motif_mask & pos_mask, "start"] = \
            match_table[motif_mask & pos_mask]["start"] + start
        match_table.loc[motif_mask & pos_mask, "end"] = \
            match_table[motif_mask & pos_mask]["start"] + length
        match_table.loc[motif_mask & pos_mask, "score_start"] = \
            match_table[motif_mask & pos_mask]["score_start"] + start
        match_table.loc[motif_mask & pos_mask, "score_end"] = \
            match_table[motif_mask & pos_mask]["score_start"] + length

        match_table.loc[motif_mask & (~pos_mask), "start"] = \
            match_table[motif_mask & (~pos_mask)]["start"] + rc_start
        match_table.loc[motif_mask & (~pos_mask), "end"] = \
            match_table[motif_mask & (~pos_mask)]["start"] + length
        match_table.loc[motif_mask & (~pos_mask), "score_start"] = \
            match_table[motif_mask & (~pos_mask)]["score_start"] + rc_start
        match_table.loc[motif_mask & (~pos_mask), "score_end"] = \
            match_table[motif_mask & (~pos_mask)]["score_start"] + length

    # For each hit, compute the total absolute importance, fraction absolute
    # importance, and IC-weighted importance average
    match_table["imp_total_score"] = np.nan
    match_table["imp_total_signed_score"] = np.nan
    match_table["imp_ic_avg_score"] = np.nan
    score_length = act_scores_matched.shape[1]
    for i, row in match_table.iterrows():
        if row["score_start"] < 0 or row["score_end"] >= score_length:
            print("Hit at %s:%d-%d is outside of importance score range" % (
                row["chrom"], row["start"], row["end"]
            ))
            continue
        scores = np.sum(act_scores_matched[
            row["example_index"], row["score_start"]:row["score_end"]
        ], axis=1)  # Flatten from L x 4 to L-array
        match_table.loc[i, "imp_total_score"] = np.sum(np.abs(scores))
        match_table.loc[i, "imp_total_signed_score"] = np.sum(scores)

        ic = ic_dict[row["key"]]
        if row["strand"] == "-":
            ic = np.flip(ic)

        match_table.loc[i, "imp_ic_avg_score"] = np.mean(ic * scores)

    # Compute the fraction importance by dividing total importance
    total_track_imp = np.sum(np.abs(act_scores_matched), axis=(1, 2))
    match_table["imp_frac_score"] = match_table["imp_total_score"] / \
        total_track_imp[match_table["example_index"]]

    # Filter out any hits that had NaN as an importance score; these hits
    # overran the importance score track boundaries
    match_table = match_table.dropna(subset=["imp_total_score"])

    # Re-order columns (and drop a few) before saving the result
    match_table = match_table[[
        "chrom", "start", "end", "key", "strand", "peak_index",
        "imp_total_signed_score", "imp_total_score", "imp_frac_score",
        "imp_ic_avg_score", "agg_sim", "mod_delta", "mod_precision",
        "mod_percentile", "fann_perclasssum_perc", "fann_perclassavg_perc"
    ]]
    match_table.to_csv(
        os.path.join(outdir, "tfm_matches.bed"), sep="\t", header=False,
        index=False
    )


if __name__ == "__main__":
    main()
