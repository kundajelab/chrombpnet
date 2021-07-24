import numpy as np
import os
import subprocess
import argparse
import h5py
import tempfile

def trim_ppm(ppm, t=0.45, min_length=3):
    maxes = np.max(ppm,-1)
    maxes = np.where(maxes>=t)

    # no base has prob>t or too small
    if len(maxes[0])==0 or maxes[0][-1]+1-maxes[0][0]<min_length:
        return None
    return ppm[maxes[0][0]:maxes[0][-1]+1]


background = np.array([0.25, 0.25, 0.25, 0.25])

def info_content(track, pseudocount=0.001):
    """
    Given an L x 4 track, computes information content for each base and
    returns it as an L-array.
    """
    num_bases = track.shape[1]
    # Normalize track to probabilities along base axis
    track_norm = (track + pseudocount) / (np.sum(track, axis=1, keepdims=True) + (num_bases * pseudocount))
    ic = track_norm * np.log2(track_norm / np.expand_dims(background, axis=0))
    return np.sum(ic, axis=1)

def write_meme_file(ppm, bg, fname):
    f = open(fname, 'w')
    f.write('MEME version 4\n\n')
    f.write('ALPHABET= ACGT\n\n')
    f.write('strands: + -\n\n')
    f.write('Background letter frequencies (from unknown source):\n')
    f.write('A %.3f C %.3f G %.3f T %.3f\n\n' % tuple(list(bg)))
    f.write('MOTIF 1 TEMP\n\n')
    f.write('letter-probability matrix: alength= 4 w= %d nsites= 1 E= 0e+0\n' % ppm.shape[0])
    for s in ppm:
        f.write('%.5f %.5f %.5f %.5f\n' % tuple(s))
    f.close()


def fetch_tomtom_matches(ppm, cwm, background=[0.25, 0.25, 0.25, 0.25], tomtom_exec_path='tomtom', motifs_db='HOCOMOCOv11_core_HUMAN_mono_meme_format.meme' ,
        n=5, trim_threshold=0.45, trim_min_length=3):
    """Fetches top matches from a motifs database using TomTom.

    Args:
        ppm: position probability matrix- numpy matrix of dimension (N,4)
        background: list with ACGT background probabilities
        tomtom_exec_path: path to TomTom executable
        motifs_db: path to motifs database in meme format
        n: number of top matches to return, ordered by p-value
        temp_dir: directory for storing temp files
        trim_threshold: the ppm is trimmed from left till first position for which
            probability for any base pair >= trim_threshold. Similarly from right.

    Returns:
        list: a list of up to n results returned by tomtom, each entry is a
            dictionary with keys 'Target ID', 'p-value', 'E-value', 'q-value'
    """

    _, fname = tempfile.mkstemp()

    # trim pfm
    #trimmed = trim_ppm(ppm, t=trim_threshold, min_length=trim_min_length)

    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * 0.3  # Cut off anything less than 20% of max score
    pass_inds = np.where(score >= trim_thresh)[0]
    trimmed = ppm[np.min(pass_inds): np.max(pass_inds) + 1]

    # can be None of no base has prob>t
    if trimmed is None:
        return []

    # trim and prepare meme file
    write_meme_file(trimmed, background, fname)

    # run tomtom
    cmd = '%s -no-ssc -oc . -verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 %s %s' % (tomtom_exec_path, fname, motifs_db)
    #print(cmd)
    out = subprocess.check_output(cmd, shell=True)

    # prepare output
    dat = [x.split('\\t') for x in str(out).split('\\n')]
    schema = dat[0]

    # meme v4 vs v5:
    if 'Target ID' in schema:
        tget_idx = schema.index('Target ID')
    else:
        tget_idx = schema.index('Target_ID')

    pval_idx, eval_idx, qval_idx =schema.index('p-value'), schema.index('E-value'), schema.index('q-value')

    r = []
    for t in dat[1:min(1+n, len(dat)-1)]:
        if t[0]=='':
            break

        mtf = {}
        mtf['Target_ID'] = t[tget_idx]
        mtf['p-value'] = float(t[pval_idx])
        mtf['E-value'] = float(t[eval_idx])
        mtf['q-value'] = float(t[qval_idx])
        r.append(mtf)

    os.system('rm ' + fname)
    return r

def fetch_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--modisco_h5py", required=True, type=str)
    parser.add_argument("-o", "--outfile", required=True, type=str)
    parser.add_argument("-d", "--meme_motif_db", required=True, type=str)
    parser.add_argument("-n", "--top_n_matches", type=int, default=3, help="Max number of matches to return from TomTom")
    parser.add_argument("-tt", "--tomtom_exec", type=str, default='tomtom')
    parser.add_argument("-th", "--trim_threshold", type=float, default=0.45, help="Trim threshold for trimming long motif, trim to those with at least prob trim_threshold on both ends")
    parser.add_argument("-tm", "--trim_min_length", type=int, default=3, help="Minimum acceptable length of motif after trimming")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args = fetch_args()
    f = h5py.File(args.modisco_h5py, 'r')

    # get pfms
    ppms = []
    cwms = []
    seqlet_tally = []

    for metacluster_name in f['metacluster_idx_to_submetacluster_results']:
        metacluster = f['metacluster_idx_to_submetacluster_results'][metacluster_name]
        if metacluster['activity_pattern'][0] == 1:
            num_patterns = len(metacluster['seqlets_to_patterns_result']['patterns'])-1

            for i in range(num_patterns):
                ppm = np.array(metacluster['seqlets_to_patterns_result']['patterns']['pattern_{}'.format(i)]['sequence']['fwd'])
                num_seqlets = len(metacluster['seqlets_to_patterns_result']['patterns']['pattern_{}'.format(i)]['seqlets_and_alnmts']['seqlets'])
                ppms.append(ppm)
                seqlet_tally.append(num_seqlets)

                cwm = np.array(metacluster['seqlets_to_patterns_result']['patterns']['pattern_{}'.format(i)]["task0_contrib_scores"]['fwd'])
                cwms.append(cwm)

    f.close()
    res = []
    for i,x in enumerate(ppms):
        res.append(fetch_tomtom_matches(x, cwms[i], tomtom_exec_path=args.tomtom_exec, motifs_db=args.meme_motif_db,
            n=args.top_n_matches, trim_threshold=args.trim_threshold, trim_min_length=args.trim_min_length))

    # write output. Skipping those patterns which disappear after trimming or have no matches
    with open(args.outfile, 'w') as f:
        # write header
        f.write("Pattern")
        f.write("\tnum_seqlets")
        for i in range(args.top_n_matches):
            f.write("\tmatch_{}\tq-value".format(i+1))
        f.write("\n")

        for i,r in enumerate(res):
            if len(r)>0:
                f.write("pattern_{}".format(i))
                f.write("\t{}".format(seqlet_tally[i]))
                for match in r:
                    f.write("\t{}\t{}".format(match['Target_ID'], match['q-value']))

                # when fewer than n matches are found
                if len(r) != args.top_n_matches:
                    f.write("\t\t"*(args.top_n_matches-len(r)))
                f.write("\n")

