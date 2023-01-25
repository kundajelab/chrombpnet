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
    with resources.path("chrombpnet.data", default_data_file_entry.value) as f:
        data_file_path=f
    return data_file_path

def print_meme_motif_file():
    with resources.path("chrombpnet.data", DefaultDataFile.motifs_meme.value) as f:
        data_file_path=f
    print(f)
    
