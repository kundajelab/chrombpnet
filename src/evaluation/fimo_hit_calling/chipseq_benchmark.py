import pandas as pd
import numpy as np
import pyBigWig as pybw

plus_strand = pybw.open("/oak/stanford/groups/akundaje/vir/tfatlas/processed_data/ENCSR257RKC/ENCSR257RKC_plus.bigWig")
minus_strand = pybw.open("/oak/stanford/groups/akundaje/vir/tfatlas/processed_data/ENCSR257RKC/ENCSR257RKC_minus.bigWig")

def get_motif_instance_profiles(motif_hit_table, chip_seq_path, profile_size=100):
    """
    Given the hit table of a specific motif, obtains the set of ChIP-nexus
    profiles of length `profile_size` centered around each motif hit's center.
    Returns an N x P x 2 array of profiles for the two strands. The first
    strand is the strand the motif is on, the second strand is the opposite.
    """
    profile_path = true_chipnexus_profile_paths[tf_name]
    
    # Get the coordinates of the profiles
    prof_coords = motif_hit_table[["chrom", "start", "end"]].values
    mids = (prof_coords[:, 1] + prof_coords[:, 2]) // 2
    prof_coords[:, 1] = mids - (profile_size // 2)
    prof_coords[:, 2] = prof_coords[:, 1] + profile_size
    
    strands = motif_hit_table["strand"].values
    
    profiles = np.empty((len(motif_hit_table), profile_size, 2))
    for i in range(len(prof_coords)):
            chrom, start, end = prof_coords[i]
            prof_plus = np.nan_to_num(plus_strand.values(chrom, start, end))
            prof_minus = np.nan_to_num(minus_strand.values(chrom, start, end))
            if strands[i] == "-":
                prof_plus = np.flip(prof_plus)
                prof_minus = np.flip(prof_minus)
            prof = np.vstack((prof_plus,prof_minus))
            profiles[i] = prof
    
    return profiles
    
def compute_footprint_scores(
    profiles, strand_0_pos, strand_1_pos, divide_by_valley=False, expand=2
):
    """
    For each profile in the N x P x 2 array `profiles`, computes a footprint
    score. For each profile, computes the sum of the profile at positions
    given by `strand_0_pos` and `strand_1_pos` (indices out of P). If
    `divide_by_valley` is provided, then this sum is divided by the position
    between the two given positions. If `expand` is nonzero, then instead of
    summing just single positions, this will sum a window by extending each
    position left and right by this amount. Returns an N-array.
    """
    assert expand >= 0
    strand_0_slice = slice(
        strand_0_pos - expand, strand_0_pos + 1 + expand 
    )
    strand_1_slice = slice(
        strand_1_pos - expand, strand_1_pos + 1 + expand 
    )
    
    strand_0_sum = np.sum(profiles[:, strand_0_slice, 0], axis=1)
    strand_1_sum = np.sum(profiles[:, strand_1_slice, 1], axis=1)
    
    if divide_by_valley:
        valley_pos = (strand_0_pos + strand_1_pos) // 2
        valley_slice = slice(
            valley_pos - expand, valley_pos + 1 + expand 
        )
        valley_0_sum = np.sum(profiles[:, valley_slice, 0], axis=1)
        valley_1_sum = np.sum(profiles[:, valley_slice, 1], axis=1)
        return (strand_0_sum / valley_0_sum) + (strand_1_sum / valley_1_sum)
    else:
        return strand_0_sum + strand_1_sum