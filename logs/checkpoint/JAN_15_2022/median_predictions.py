import h5py
import numpy as np
import matplotlib.pyplot as plt

main="results/chrombpnet/ATAC_PE/GM12878/ATAC_PE_12.30.2021/chrombpnet_model/"
cwbm = main+"chrombpnet_wo_bias_predictions.h5"
cm = main + "chrombpnet_predictions.h5"
bm = main + "bias_predictions.h5"

print("peaks and non peaks")


bm_data = h5py.File(bm,"r")
print(len(bm_data["predictions"]["logcounts"]))
mse = np.median(bm_data["predictions"]["logcounts"])
print("bias mse", mse)

cwbm_data = h5py.File(cwbm,"r")
mse = np.median(cwbm_data["predictions"]["logcounts"])
print("chrombpnet wo bias mse", mse)

cm_data = h5py.File(cm,"r")
mse = np.median(cm_data["predictions"]["logcounts"])
print("chrombpnet mse", mse)

print("only peaks")


bm_data = h5py.File(bm,"r")
print(np.array(bm_data["predictions"]["logcounts"])[np.array(bm_data["coords"]["coords_peak"])==1].shape)
mse = np.median(bm_data["predictions"]["logcounts"][np.array(bm_data["coords"]["coords_peak"])==1])
print("bias mse", mse)

cwbm_data = h5py.File(cwbm,"r")
mse = np.median(cwbm_data["predictions"]["logcounts"][np.array(cwbm_data["coords"]["coords_peak"])==1])
print("chrombpnet wo bias mse", mse)

cm_data = h5py.File(cm,"r")
mse = np.median(cm_data["predictions"]["logcounts"][np.array(cm_data["coords"]["coords_peak"])==1])
print("chrombpnet mse", mse)


print("only non peaks")


bm_data = h5py.File(bm,"r")
print(np.array(bm_data["predictions"]["logcounts"])[np.array(bm_data["coords"]["coords_peak"])==0].shape)
mse = np.median(bm_data["predictions"]["logcounts"][np.array(bm_data["coords"]["coords_peak"])==0])
print("bias mse", mse)

cwbm_data = h5py.File(cwbm,"r")
mse = np.median(cwbm_data["predictions"]["logcounts"][np.array(cwbm_data["coords"]["coords_peak"])==0])
print("chrombpnet wo bias mse", mse)

cm_data = h5py.File(cm,"r")
mse = np.median(cm_data["predictions"]["logcounts"][np.array(cm_data["coords"]["coords_peak"])==0])
print("chrombpnet mse", mse)


plt.scatter(bm_data["predictions"]["logcounts"][np.array(bm_data["coords"]["coords_peak"])==0], cwbm_data["predictions"]["logcounts"][np.array(cwbm_data["coords"]["coords_peak"])==0])
plt.show()
