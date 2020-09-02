import numpy as np

def get_vals_from_gkm_line(seq):
    return np.asarray([[float(i) for i in i.split(',')] for i in seq.split(';')])
