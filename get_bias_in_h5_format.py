import tensorflow as tf
import pandas as pd

data = pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_atac.csv", sep=',', header=None)
data = data[data[1]=="HEPG2"].reset_index(drop=True)

print(data.head())

for i,r in data.iterrows():
	ppath = r[2]+"/bias_model/new_model_formats/bias/"
	print(ppath)
	modelf = tf.keras.models.load_model(ppath)
	odir="/oak/stanford/groups/akundaje/anusri/dil_layer_expts/HEPG2/"+r[2].strip().split("/")[-1]
	print(odir)
	modelf.save(odir+"/bias_model/bias.h5")
