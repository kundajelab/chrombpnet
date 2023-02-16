import os
import subprocess
import numpy as np
import pandas as pd
import h5py
import utils
import pyBigWig
import pyfaidx
from tqdm import tqdm
from scipy import signal

def dna_to_one_hot(seqs):
    """
    Converts a list of DNA ("ACGT") sequences to one-hot encodings, where the
    position of 1s is ordered alphabetically by "ACGT". `seqs` must be a list
    of N strings, where every string is the same length L. Returns an N x L x 4
    NumPy array of one-hot encodings, in the same order as the input sequences.
    All bases will be converted to upper-case prior to performing the encoding.
    Any bases that are not "ACGT" will be given an encoding of all 0s.
    """
    seq_len = len(seqs[0])
    assert np.all(np.array([len(s) for s in seqs]) == seq_len)

    # Join all sequences together into one long string, all uppercase
    seq_concat = "".join(seqs).upper() + "ACGT"
    # Add one example of each base, so np.unique doesn't miss indices later

    one_hot_map = np.identity(5)[:, :-1].astype(np.int8)

    # Convert string into array of ASCII character codes;
    base_vals = np.frombuffer(bytearray(seq_concat, "utf8"), dtype=np.int8)

    # Anything that's not an A, C, G, or T gets assigned a higher code
    base_vals[~np.isin(base_vals, np.array([65, 67, 71, 84]))] = 85

    # Convert the codes into indices in [0, 4], in ascending order by code
    _, base_inds = np.unique(base_vals, return_inverse=True)

    # Get the one-hot encoding for those indices, and reshape back to separate
    return one_hot_map[base_inds[:-4]].reshape((len(seqs), seq_len, 4))


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

def export_motifs_meme_format(pfms, out_dir):
    """
    Exports motifs to an output directory as PFMs for fimo.
    Arguments:
        `pfms`: a dictionary mapping keys to N x 4 NumPy arrays (N may be
            different for each PFM); `{key}.pfm` will be the name of each saved
            motif
        `out_dir`: directory to save motifs in a meme format
    """
    outfile = os.path.join(out_dir, "motifs.pfm")
    with open(outfile, "w") as f:
        f.write("MEME version 4" + "\n"+ "\n")
        f.write("ALPHABET= ACGT" + "\n"+ "\n")
        f.write("Background letter frequencies \n"+"A 0.25 C 0.25 G 0.25 T 0.25 \n"+ "\n")
    f.close()

    for key, pfm in pfms.items():
        with open(outfile, "a") as f:
            f.write(f"MOTIF {key}" + "\n")
            f.write("letter-probability matrix:"+ "\n")
            for i in range(pfm.shape[0]):
                f.write(" ".join([str(x) for x in pfm[i,:]]) + "\n")
            f.write("\n")
    f.close()

def run_moods(out_dir, reference_fasta, peaks_bed_path, pval_thresh=0.01):
    """
    Runs MOODS on every `.pfm` file in `out_dir`. Outputs the results for each
    PFM into `out_dir/moods_out.csv`.
    Arguments:
        `out_dir`: directory with PFMs
        `reference_fasta`: path to reference Fasta to use
        `pval_thresh`: threshold p-value for MOODS to use
    """ 
    comm = ["bedtools"]
    comm += ["getfasta"]
    comm += ["-fi"]
    comm += [reference_fasta]
    comm += ["-bed"]
    comm += [peaks_bed_path]
    comm += ["-fo"]
    comm += [os.path.join(out_dir, "peaks.fa")]
    proc = subprocess.Popen(comm)
    proc.wait()
    
    pfm_files = [p for p in os.listdir(out_dir) if p.endswith(".pfm")]
    comm = ["moods-dna.py"]
    comm += ["--batch"]
    comm += ["-m"]
    comm += [os.path.join(out_dir, pfm_file) for pfm_file in pfm_files]
    comm += ["-s", os.path.join(out_dir, "peaks.fa")]
    comm += ["-p", str(pval_thresh)]
    comm += ["-o", os.path.join(out_dir, "moods_out.csv")]
    proc = subprocess.Popen(comm)
    proc.wait()
    
def run_fimo(out_dir, reference_fasta, peaks_bed_path, pval_thresh=0.01):


	comm = ["bedtools"]
	comm += ["getfasta"]
	comm += ["-fi"]
	comm += [reference_fasta]
	comm += ["-bed"]
	comm += [peaks_bed_path]
	comm += ["-fo"]
	comm += [os.path.join(out_dir, "peaks.fa")]
	proc = subprocess.Popen(comm)
	proc.wait()
    
	comm = ["fimo"]
	comm += ["--thresh", str(pval_thresh)]
	comm += ["--no-qvalue", "--verbosity", "5"]
	comm += ["--max-stored-scores", "1000000"]
	comm += ["--parse-genomic-coord"]
	comm += ["--bfile", "--uniform--"]
	comm += ["--o", os.path.join(out_dir, "fimo_out")]
	comm += [os.path.join(out_dir, "motifs.pfm")]
	comm += [os.path.join(out_dir, "peaks.fa")]
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
            hit_pos=int(tokens[0].split(":")[1].split("-")[0])+int(tokens[2])
            g.write("\t".join([
                tokens[0].split(":")[0], str(hit_pos),
                str(hit_pos + len(tokens[5])), tokens[1][:-4], tokens[3],
                tokens[4]
            ]) + "\n")
#             g.write("\t".join([
#                 tokens[0].split()[0], tokens[2],
#                 str(int(tokens[2]) + len(tokens[5])), tokens[1][:-4], tokens[3],
#                 tokens[4]
#             ]) + "\n")
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

def fimo_hits_to_bed(moods_out_csv_path, moods_out_bed_path):
    """
    Converts MOODS hits into BED file.
    """
    f = pd.read_csv(moods_out_csv_path, sep="\t", header=0)
    print(f.shape)
    
    cols = ["sequence_name", "start", "stop", "motif_id", "strand", "score" ]
    g = f[cols]
    print(g.head())
    g = g[g['start'].notna()]
    g = g[g['stop'].notna()]
    print(g.shape)
    g["start"] = g["start"].astype(int)
    g["stop"] = g["stop"].astype(int)

    g.to_csv(moods_out_bed_path, sep="\t", header=False, index=False)



def filter_hits_for_peaks(
    moods_out_bed_path, filtered_hits_path, peak_bed_path, num_cols_in_peaks
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
        #print(tokens)
        g.write("\t".join((tokens[:-num_cols_in_peaks])))
        if num_cols_in_peaks==3:
        	peak_index = peak_index_map["%s:%s-%s" % tuple(tokens[-num_cols_in_peaks:])]
        if num_cols_in_peaks==10:
        	peak_index = peak_index_map["%s:%s-%s" % tuple(tokens[-num_cols_in_peaks:-num_cols_in_peaks+3])]
        
        g.write("\t" + peak_index + "\n")
    f.close()
    g.close()


def compute_hits_importance_scores(
    hits_bed_path, shap_scores_hdf5_path, peak_table, peak_bed_path, out_path, method
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

    imp_scores = utils.import_shap_scores_from_bigwig(
        shap_scores_hdf5_path, peak_table
    )
    
    peak_table = pd.read_csv(
        peak_bed_path, sep="\t", header=None, index_col=False,
        usecols=[0, 1, 2], names=["peak_chrom", "peak_start", "peak_end"]
    )

    hit_table = pd.read_csv(
        hits_bed_path, sep="\t", header=None, index_col=False,
        names=["chrom", "start", "end", "key", "strand", "score", "peak_index"]
    )

    
    # Merge in the peak starts/ends to the hit table
    merged_hits = pd.merge(
        hit_table, peak_table, left_on="peak_index", right_index=True
    )

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
            if method=="sum_norm":
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
            if method=="mean":
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


def get_hits(
    pfm_dict, reference_fasta, peak_bed_path, shap_scores_hdf5_path, peak_table_in,
    moods_pval_thresh=0.01, temp_dir=None, method="sum", hit_caller="fimo"
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

    pfm_keys = list(pfm_dict.keys())
    
    if hit_caller == "fimo":
        export_motifs_meme_format(pfm_dict, temp_dir)
        print("Running fimo...")
        run_fimo(temp_dir, reference_fasta, peak_bed_path, pval_thresh=moods_pval_thresh)
        fimo_hits_to_bed(
        	os.path.join(temp_dir, "fimo_out/fimo.tsv"),
        	os.path.join(temp_dir, "moods_out.bed")
        )
    
    if hit_caller == "moods":
    	export_motifs(pfm_dict, temp_dir)
    	print("Running moods...")
    	run_moods(temp_dir, reference_fasta, peak_bed_path, pval_thresh=moods_pval_thresh)
    	
    	(
    		os.path.join(temp_dir, "moods_out.csv"),
    		os.path.join(temp_dir, "moods_out.bed")
    	)
    

    print("Compute importance scores in hits")
    filter_hits_for_peaks(
        os.path.join(temp_dir, "moods_out.bed"),
        os.path.join(temp_dir, "moods_peaks_annotate.bed"),
        peak_bed_path, 3
    )
    
    compute_hits_importance_scores(
        os.path.join(temp_dir, "moods_peaks_annotate.bed"),
        shap_scores_hdf5_path, peak_table_in, peak_bed_path,
        os.path.join(temp_dir, "moods_scored.bed"), method
    )

    hit_table = import_moods_hits(
        os.path.join(temp_dir, "moods_scored.bed")
    )

    return hit_table

def resolve_overlaps(output_dir,filtered_hits_path,contribs_bw,reference_fasta, pfms):

	comm = ["bedtools", "sort", "-i"]
	comm += [filtered_hits_path]
	comm += ["|"]
	comm += ["bedtools", "cluster", "-i", "stdin"]
	f = open(os.path.join(output_dir, "moods_filtered_clustered.bed"), "w")
	proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
	proc.wait()
	f.close()
	
	data = pd.read_csv(os.path.join(output_dir, "moods_filtered_clustered.bed"), sep="\t", header=None)
	bw = pyBigWig.open(contribs_bw)
	genome = pyfaidx.Fasta(reference_fasta)
	
	corrs=[]
	for i,r in tqdm(data.iterrows()):
		chrom = r[0]
		start = r[1]
		end = r[2]
		key = r[3]
		
		val = np.nan_to_num(bw.values(chrom,start,end)).reshape((-1,1))
		seq = dna_to_one_hot([str(genome[chrom][start:end])])[0]
		output = val*seq
		
		# check correlation with both orginal pfm motif and is reverse complement
		corr1 = np.max(signal.correlate2d(np.abs(pfms[key][::-1, ::-1]),output, mode="valid"))
		corr2 = np.max(signal.correlate2d(np.abs(pfms[key]),output, mode="valid"))

		corr = np.max([corr1, corr2])
		corrs.append(corr/(end-start)) # dividing it by length - this prefers shorter motifs while resolving clashes - so if you have 2 CTCF variants in your list and one is shorter that will get most of the hits
		
	data[9] = corrs
			
	new_data = []
	clusters = list(set(data[8].values))
	for idx in tqdm(clusters):
		obc = data[data[8]==idx].reset_index(drop=True)
		if len(obc) == 1:
			new_data.append(obc.loc[0,:].values)
		elif len(obc) == 0:
			continue
		else:
			while len(obc) > 0:
				idx = obc[9].idxmax()
				new_data.append(obc.loc[idx,:].values)
				width = (obc.loc[idx,2] - obc.loc[idx,1])//4
				obc = obc[ (obc[2] < (obc.loc[idx,1]+width)) | (obc[1] > (obc.loc[idx,2]-width))].reset_index(drop=True)
		
	new_data = pd.DataFrame(new_data)

	return new_data