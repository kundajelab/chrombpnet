import tensorflow as tf
from tensorflow.keras.models import load_model
import argparse
from tensorflow.keras.utils import get_custom_objects
import os
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Add, Concatenate, Lambda, Flatten

parser = argparse.ArgumentParser(description="converting model types")
parser.add_argument("-i", "--input_model")
parser.add_argument("-o","--output_dir")
parser.add_argument("-f","--file_path")

args = parser.parse_args()


output_path=os.path.join(args.output_dir, args.file_path+"/")
#custom_objects={"tf": tf}
custom_objects={"logcount_predictions":  Lambda(lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)), "tf": tf}
get_custom_objects().update(custom_objects)
model=load_model(args.input_model,compile=False)
model.save(output_path)


command = "cd "+args.output_dir+" && tar -cf "+args.file_path+".tar "+args.file_path+"/"
#command="tar -cf "+args.output_path+".tar"+" "+args.output_path
os.system(command)


#tf.saved_model.save(model,args.output_path)
