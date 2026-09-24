import numpy as np
from chrombpnet.training.utils.data_utils import one_hot


def get_seq(peaks_df, genome, width):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []
    peaks_used = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        if len(sequence) == width:
            vals.append(sequence)
            peaks_used.append(True)
        else:
            peaks_used.append(False)

    return one_hot.dna_to_one_hot(vals), np.array(peaks_used)


def load_model_wrapper(args):
    # read .h5 model (args.model_h5, or a path)
    from chrombpnet.training.utils import model_io
    model_h5 = args if isinstance(args, str) else args.model_h5
    model = model_io.load_model(model_h5, compile=False)
    print("got the model")
    model.summary()
    return model
