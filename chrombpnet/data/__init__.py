from importlib import resources
from enum import Enum

class DefaultDataFile(Enum):
    atac_ref_motifs = "ATAC.ref.motifs.txt"
    dnase_ref_motifs = "DNASE.ref.motifs.txt"
    motif_to_pwm_atac = "motif_to_pwm.ATAC.tsv"
    motif_to_pwm_dnase = "motif_to_pwm.DNASE.tsv"
    motif_to_pwm_tf = "motif_to_pwm.TF.tsv"
    motifs_meme = "motifs.meme.txt"
    

def get_default_data_path(default_data_file_entry):    
    # a pathlib.Path: the package is installed as regular files (no zip imports)
    data_file_path = resources.files("chrombpnet.data") / default_data_file_entry.value
    return data_file_path

def print_meme_motif_file():
    data_file_path = get_default_data_path(DefaultDataFile.motifs_meme)
    print(data_file_path)
