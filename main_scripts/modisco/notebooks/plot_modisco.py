import numpy as np
import modisco
import h5py
from matplotlib import pyplot as plt
import modisco.visualization
from modisco.visualization import viz_sequence
import seaborn as sns

def plot_modisco(path):
    hdf5_results = h5py.File(path,'r')

    print("Metaclusters heatmap")
    
    activity_patterns = np.array(hdf5_results['metaclustering_results']['attribute_vectors'])[
                        np.array(
            [x[0] for x in sorted(
                    enumerate(hdf5_results['metaclustering_results']['metacluster_indices']),
                   key=lambda x: x[1])])]
    sns.heatmap(activity_patterns, center=0)
    plt.show()

    metacluster_names = [
        x.decode("utf-8") for x in 
        list(hdf5_results["metaclustering_results"]
             ["all_metacluster_names"][:])]

    all_patterns = []

    for metacluster_name in metacluster_names:
        print(metacluster_name)
        metacluster_grp = (hdf5_results["metacluster_idx_to_submetacluster_results"]
                                       [metacluster_name])
        if metacluster_grp["activity_pattern"][0] == 1:
            print("activity pattern:",metacluster_grp["activity_pattern"][:])
            all_pattern_names = [x.decode("utf-8") for x in 
                                    list(metacluster_grp["seqlets_to_patterns_result"]
                                                 ["patterns"]["all_pattern_names"][:])]
            if (len(all_pattern_names)==0):
                print("No motifs found for this activity pattern")
            for pattern_name in all_pattern_names:
                print(metacluster_name, pattern_name)
                all_patterns.append((metacluster_name, pattern_name))
                pattern = metacluster_grp["seqlets_to_patterns_result"]["patterns"][pattern_name]
                print("total seqlets:",len(pattern["seqlets_and_alnmts"]["seqlets"]))
                background = np.array([0.27, 0.23, 0.23, 0.27])
                print("Hypothetical scores:")
                viz_sequence.plot_weights(pattern["task0_hypothetical_contribs"]["fwd"])
                print("Actual importance scores:")
                viz_sequence.plot_weights(pattern["task0_contrib_scores"]["fwd"])
                print("onehot, fwd and rev:")
                viz_sequence.plot_weights(viz_sequence.ic_scale(np.array(pattern["sequence"]["fwd"]),
                                                                background=background)) 
                viz_sequence.plot_weights(viz_sequence.ic_scale(np.array(pattern["sequence"]["rev"]),
                                                                background=background))

    hdf5_results.close()
