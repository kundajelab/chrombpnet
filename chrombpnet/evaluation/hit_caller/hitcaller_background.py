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
    os.makedirs(os.path.join(out_dir, "motifs/"), exist_ok=True) 
    
    out_dir_new = os.path.join(out_dir, "motifs/")
    for key, pfm in pfms.items():
    	outfile = os.path.join(out_dir_new, key+".pfm")
    	with open(outfile, "w") as f:
        	f.write("MEME version 4" + "\n"+ "\n")
        	f.write("ALPHABET= ACGT" + "\n"+ "\n")
        	f.write("Background letter frequencies \n"+"A 0.25 C 0.25 G 0.25 T 0.25 \n"+ "\n")
    	f.close()

    	with open(outfile, "a") as f:
            f.write(f"MOTIF {key}" + "\n")
            f.write("letter-probability matrix:"+ "\n")
            for i in range(pfm.shape[0]):
                f.write(" ".join([str(x) for x in pfm[i,:]]) + "\n")
            f.write("\n")
    	f.close()

    
def run_fimo(out_dir, pval_thresh=0.01):

	os.makedirs(os.path.join(out_dir, "fimo_out/"), exist_ok=True) 

	pvals_per_motif = {}
	for file in os.listdir(os.path.join(out_dir,"motifs/")):
		
		found_enough_hits=False
		motif_pval = float(pval_thresh)
		if file.endswith(".pfm"):
			while not found_enough_hits:
				print(os.path.join(os.path.join(out_dir,"motifs/"), file))
				file_path = os.path.join(os.path.join(out_dir,"motifs/"), file)
				comm = ["fimo"]
				comm += ["--thresh", str(motif_pval)]
				comm += ["--no-qvalue", "--verbosity", "5"]
				comm += ["--max-stored-scores", "100000000"]
				comm += ["--bfile", os.path.join(out_dir, "frequency.txt")]
				comm += ["--o", os.path.join(out_dir, "fimo_out/"+file.replace(".pfm",""))]
				comm += [file_path]
				comm += [os.path.join(out_dir, "peaks.fa")]
				proc = subprocess.Popen(comm)
				proc.wait()
		
				if motif_pval == 1:
					pvals_per_motif[file] = 1
					print(file, motif_pval)
					break
				
				print(os.path.join(out_dir, "fimo_out/"+file.replace(".pfm","")+"/fimo.tsv"))		
				if  os.path.join(out_dir, "fimo_out/"+file.replace(".pfm","")+"/fimo.tsv"):

					f = pd.read_csv(os.path.join(out_dir, "fimo_out/"+file.replace(".pfm","")+"/fimo.tsv"), sep="\t", header=0)
					print(f.shape)
					cols = ["sequence_name", "start", "stop", "motif_id", "strand", "score" ]
					try:
						g = f[cols]
						print(g.shape)
						if g.shape[0] >= 700000:
							found_enough_hits=True
							pvals_per_motif[file] = motif_pval
							print(file, motif_pval)
						else:
							motif_pval = motif_pval*10
							os.system("rm -r "+os.path.join(out_dir, "fimo_out/"+file.replace(".pfm","")))
					except:
						motif_pval = motif_pval*10
						os.system("rm -r "+os.path.join(out_dir, "fimo_out/"+file.replace(".pfm","")))

					if motif_pval > 1:
						motif_pval = 1
					
	return pvals_per_motif				


def fimo_hits_to_bed(moods_out_path, moods_out_bed_path):
	"""
	Converts MOODS hits into BED file.
	"""
	frames = []
	for file in os.listdir(moods_out_path):
	
		#if file.endswith(".tsv"):
		#	print(file)
		odir = os.path.join(moods_out_path,file)
		print(odir)
		if os.path.isfile(os.path.join(odir,"fimo.tsv")):
			print(os.path.join(odir,"fimo.tsv"))		

			f = pd.read_csv(os.path.join(odir,"fimo.tsv"), sep="\t", header=0)

			cols = ["sequence_name", "start", "stop", "motif_id", "strand", "score" ]
			g = f[cols]
			#print(g.head())
			g = g[g['start'].notna()]
			g = g[g['stop'].notna()]
			#print(g.shape)
			g["start"] = g["start"].astype(int)
			g["stop"] = g["stop"].astype(int)

			frames.append(g)
    	
	new_data = pd.concat(frames)
	new_data.to_csv(moods_out_bed_path, sep="\t", header=False, index=False)



def compute_hits_importance_scores(
    hits_bed_path, imp_scores, out_path, method
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

    hit_table = pd.read_csv(
        hits_bed_path, sep="\t", header=None, index_col=False,
        names=["peak_index", "start", "stop", "key", "strand", "score"]
    )
    
    hit_table["peak_index"] = hit_table["peak_index"].astype(int)
    # Compute start and end of each motif relative to the peak
    hit_table["motif_rel_start"] = \
        hit_table["start"] - 1 # 1 - coordinate system
    hit_table["motif_rel_end"] = \
         hit_table["stop"] # 1 - coordinate system

    scores = np.empty(len(hit_table))
    for peak_index, group in hit_table.groupby("peak_index"):
        # Iterate over grouped table by peak
        score_track = np.abs(imp_scores[peak_index])
        total_score = np.abs(np.mean(score_track))
        for i, row in group.iterrows():
            if method=="sum_norm":
                scores[i] = np.abs(np.sum(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                ))/total_score
            if method=="sum":
                scores[i] = np.abs(np.sum(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                ))
            if method=="mean_norm":
                scores[i] = np.abs(np.mean(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                ))/total_score
                #scores[i] = np.mean(
                #        score_track[row["motif_rel_start"]:row["motif_rel_end"]]/total_score
                #)
            if method=="mean":
                scores[i] = np.abs(np.mean(
                        score_track[row["motif_rel_start"]:row["motif_rel_end"]]
                ))
                
    hit_table["imp_frac_score"] = scores
    new_hit_table = hit_table[[
        "peak_index", "start", "stop", "key", "strand", "score", "imp_frac_score"
    ]]
    new_hit_table.to_csv(out_path, sep="\t", header=False, index=False)
    return new_hit_table


def get_hits(
    pfm_dict, shap_scores_hdf5_path, moods_pval_thresh=0.01, temp_dir=None, method="sum", hit_caller="fimo"
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
    
    if not os.path.exists(temp_dir):
    	os.makedirs(temp_dir)

    pfm_keys = list(pfm_dict.keys())
    
    import deepdish
    
    scores = deepdish.io.load(shap_scores_hdf5_path)
    
    num_seqs = scores["projected_shap"]["seq"].shape[0]
    imp_scores_array = np.sum(scores["projected_shap"]["seq"],axis=1)
    print(imp_scores_array.shape)
    
    bases = np.array(["A", "C", "G", "T", "N"])
    #one_hot = np.where(scores["raw"]["seq"]==1,1,0)
    #one_hot = scores["raw"]["seq"]
    one_hot = scores["raw"]["seq"]
    one_hot = np.transpose(one_hot,(0,2,1))
    batch_inds, seq_inds, base_inds = np.where(one_hot)
    
    one_hot_inds = np.tile(one_hot.shape[2], one_hot.shape[:2])
    one_hot_inds[batch_inds, seq_inds] = base_inds
    seq_array = bases[one_hot_inds]
    sequences = ["".join(seq) for seq in seq_array]
    
    assert(len(sequences) == imp_scores_array.shape[0])
    ofile=open(os.path.join(temp_dir,"peaks.fa"),"w")
    for idx in range(len(sequences)):
    	ofile.write(">"+str(idx))
    	ofile.write("\n")
    	ofile.write(sequences[idx])
    	ofile.write("\n")
    ofile.close()
    
    comm = ["fasta-get-markov"]
    comm += [os.path.join(temp_dir, "peaks.fa")]
    comm += [os.path.join(temp_dir, "frequency.txt")]
    comm += ["-m", "1"]
    proc = subprocess.Popen(comm)
    proc.wait()	
	 	   
    if hit_caller == "fimo":
        export_motifs_meme_format(pfm_dict, temp_dir)
        print("Running fimo...")
        pvals_per_motif = run_fimo(temp_dir, pval_thresh=moods_pval_thresh)
        fimo_hits_to_bed(
        	os.path.join(temp_dir, "fimo_out/"),
        	os.path.join(temp_dir, "moods_out.bed")
        )
    
    print("Compute importance scores in hits")
    
    hit_table = compute_hits_importance_scores(
        os.path.join(temp_dir, "moods_out.bed"),
        imp_scores_array,
        os.path.join(temp_dir, "moods_scored.bed"), method
    )


    return hit_table, pvals_per_motif

