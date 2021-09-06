import matplotlib
matplotlib.use('pdf')
import h5py
import numpy as np
import modisco
import sys
import pickle

infile = sys.argv[1]
outfile = sys.argv[2]
score_type = sys.argv[3]
cores = int(sys.argv[4])
seqlets = int(sys.argv[5])

open_file = open(infile, 'rb')
score_dict = pickle.load(open_file)

hyp_impscores = []
impscores = []
onehot_seqs = []

for i in score_dict[score_type]:
    hyp_impscores.append(score_dict[score_type][i])
    impscores.append(score_dict[score_type][i])
    onehot_seqs.append(score_dict['seq'][i])

hyp_impscores = np.array(hyp_impscores)
onehot_seqs = np.array(onehot_seqs)
impscores = hyp_impscores * onehot_seqs
null_per_pos_scores = modisco.coordproducers.LaplaceNullDist(num_to_samp=5000)

tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                    max_seqlets_per_metacluster=seqlets,
                    sliding_window_size=21,
                    flank_size=10,
                    target_seqlet_fdr=0.05,
                    min_passing_windows_frac=0.03,
                    seqlets_to_patterns_factory=modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                        n_cores=cores,
                        trim_to_window_size=30,
                        initial_flank_to_add=10,
                        final_min_cluster_size=30))(task_names=["task0"],
                            contrib_scores={'task0': impscores},
                            hypothetical_contribs={'task0': hyp_impscores},
                            null_per_pos_scores=null_per_pos_scores,
                            one_hot=onehot_seqs)

h5f = h5py.File(outfile, 'w')
tfmodisco_results.save_hdf5(h5f)
h5f.close()

