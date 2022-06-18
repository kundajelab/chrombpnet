import h5py
import pyBigWig
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import numpy as np 
from tqdm import tqdm
import scipy.ndimage

file1= "/oak/stanford/groups/akundaje/arpitas/thyroid_normal/croo/chrombpnet/models/chrombpnet_default_bias_model/chrombpnet_predictions.h5"
bigwig="/oak/stanford/groups/akundaje/arpitas/thyroid_normal/croo/chrombpnet/data/downloads/merged-bigwig.bw"
rep1=pyBigWig.open("/mnt/lab_data2/anusri/thyroid/rep1.bw")
rep2=pyBigWig.open("/mnt/lab_data2/anusri/thyroid/rep2.bw")
data = h5py.File(file1)

profile_predictions = data["predictions"]["profs"].value
prof_coords_chr = data["coords"]["coords_chrom"].value
prof_coords_center = data["coords"]["coords_center"].value
prof_coords_peak = data["coords"]["coords_peak"].value
num_bins=100
plt.rcParams["figure.figsize"]=8,8


#print(prof_coords_chr.shape)
jsd_pw=[]
jsd_rnd=[]
jsd_rep=[]

bw = pyBigWig.open(bigwig) 
for idx in tqdm(range(prof_coords_chr.shape[0])):
	chr = prof_coords_chr[idx]
	start = prof_coords_center[idx] - 500
	end = prof_coords_center[idx] + 500
	if prof_coords_peak[idx] == 0:
		continue

	pseudocount=0.001

	true_counts = np.nan_to_num(bw.values(chr,start,end ))+pseudocount
	shuffled_labels=np.random.permutation(true_counts)
	pred_probs = profile_predictions[idx]

	rep1_counts = np.nan_to_num(rep1.values(chr,start,end ))+pseudocount
	rep2_counts = np.nan_to_num(rep2.values(chr,start,end ))+pseudocount

	#true_counts = scipy.ndimage.gaussian_filter1d(true_counts, 7,axis=0, truncate=(80 / 14))
	#pred_probs = scipy.ndimage.gaussian_filter1d(pred_probs, 7,axis=0, truncate=(80 / 14))

	cur_jsd=jensenshannon(true_counts/(np.nansum(true_counts)),pred_probs)
	jsd_pw.append(cur_jsd)



	#shuffled_labels = scipy.ndimage.gaussian_filter1d(shuffled_labels, 7,axis=0, truncate=(80 / 14))
	shuffled_labels_prob=shuffled_labels/(np.nansum(shuffled_labels))


	#rep1_counts = scipy.ndimage.gaussian_filter1d(rep1_counts, 7,axis=0, truncate=(80 / 14))
	#rep2_counts = scipy.ndimage.gaussian_filter1d(rep2_counts, 7,axis=0, truncate=(80 / 14))
	curr_jsd_rnd=jensenshannon(true_counts/(np.nansum(true_counts)),shuffled_labels_prob)
	jsd_rnd.append(curr_jsd_rnd)

	rep_jsd=jensenshannon(rep1_counts/(np.nansum(rep1_counts)),rep2_counts/(np.nansum(rep2_counts)))
	jsd_rep.append(rep_jsd)




print(np.nanmedian(jsd_pw))
print(np.nanmedian(jsd_rnd))
print(np.nanmedian(jsd_rep))

#print(jsd_pw)
plt.figure()
n,bins,patches=plt.hist(jsd_pw,num_bins,facecolor='blue',alpha=0.5,label="Predicted vs Labels")
n1,bins1,patches1=plt.hist(jsd_rnd,num_bins,facecolor='black',alpha=0.5,label='Shuffled Labels vs Labels')
n1,bins1,patches1=plt.hist(jsd_rep,num_bins,facecolor='red',alpha=0.5,label='Pseudoreplicate performance')
plt.xlabel('Jensen Shannon Distance Profile Labels and Predictions in Probability Space')
plt.title("JSD Dist: ")
plt.legend(loc='best')
plt.savefig("temp.no.smooth.jsd.png",format='png',dpi=300)
