import pdb 
import argparse
import pickle
import tensorflow
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import kerasAC 
from kerasAC.generators.tiledb_predict_generator import *
from kerasAC.tiledb_config import *
from kerasAC.interpret.deepshap import * 
from kerasAC.interpret.profile_shap import * 
from kerasAC.helpers.transform_bpnet_io import *
#load the model!
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import * 
from kerasAC.custom_losses import * 
from tensorflow.keras.models import load_model, model_from_json
#from bpnet_shap_weighted import *
from keras_genomics.layers.convolutional import RevCompConv1D
import os
import numpy as np
import random as rn
import tensorflow as tf
import keras

os.environ['PYTHONHASHSEED'] = '0'
seed=1234
np.random.seed(seed)
tf.random.set_seed(seed)
rn.seed(seed)

def load_model_weights(weight_file,model):
    if weight_file is None:
        #nothing to do
        return model
    #sync the model locally if it's on AWS
    if weight_file.startswith("s3://"):
        s3_model_weights=download_s3_file(weight_file)
    else:
        s3_model_weights=weight_file
    import  h5py
    try:
        model.load_weights(s3_model_weights)#,by_name=True)
    except:
        with h5py.File(s3_model_weights,'r') as file:
            weight_file=file['model_1']
            for layer in model.layers:
                try:
                    layer_weights=weight_file[layer.name]
                except:
                    print('no weights files saved for layer:'+str(layer.name))
                    continue
                try:
                    weights = []
                    # Extract weights
                    for term in layer_weights:
                        if isinstance(layer_weights[term], h5py.Dataset):
                            # Convert weights to numpy array and prepend to list
                            weights.insert(0, np.array(layer_weights[term]))
                    # Load weights to model
                    layer.set_weights(weights)
                    print("loaded weights for layer:"+str(layer.name))
                except Exception as e:
                    print("Error: Could not load weights for layer:"+str(layer.name))
                    raise e
    return model



def parse_args():
    parser=argparse.ArgumentParser(description="wrapper to make it easier to deepSHAP bpnet models")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--chipseq",action='store_true',default=False) 
    parser.add_argument("--tfchipseq",action='store_true',default=False) 
    parser.add_argument("--model_hdf5",default=None)
    parser.add_argument("--bed_regions")
    parser.add_argument("--bed_regions_center",choices=['summit','center'])
    parser.add_argument("--tdb_array")
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--tdb_output_datasets",nargs="+")
    parser.add_argument("--tdb_input_datasets",nargs="+")
    parser.add_argument("--batch_size",type=int)
    parser.add_argument("--tdb_output_source_attribute",nargs="+",help="tiledb attribute for use in label generation i.e. fc_bigwig")
    parser.add_argument("--tdb_output_flank",nargs="+",help="flank around bin center to use in generating outputs")
    parser.add_argument("--tdb_output_aggregation",nargs="+",help="method for output aggregation; one of None, 'avg','max'")
    parser.add_argument("--tdb_output_transformation",nargs="+",help="method for output transformation; one of None, 'log','log10','asinh'")
    parser.add_argument("--tdb_input_source_attribute",nargs="+",help="attribute to use for generating model input, or 'seq' for one-hot-encoded sequence")
    parser.add_argument("--tdb_input_flank",nargs="+",help="length of sequence around bin center to use for input")
    parser.add_argument("--tdb_input_aggregation",nargs="+",help="method for input aggregation; one of 'None','avg','max'")
    parser.add_argument("--tdb_input_transformation",nargs="+",help="method for input transformation; one of None, 'log','log10','asinh'")
    parser.add_argument("--num_inputs",type=int)
    parser.add_argument("--num_outputs",type=int) 
    parser.add_argument("--out_pickle")
    parser.add_argument("--json_string")
    parser.add_argument("--weights")
    parser.add_argument("--task_index",type=int)
    parser.add_argument("--num_threads",type=int)
    return parser.parse_args() 


class BiasLayer(keras.layers.Layer):
    def __init__(self, *args, **kwargs):
        super(BiasLayer, self).__init__(*args, **kwargs)

    def build(self, input_shape):
        self.bias = self.add_weight('bias',
                                    shape=input_shape[1:],
                                    initializer='zeros',
                                    trainable=True)
    def call(self, x):
        return x + self.bias
	
def load_model_wrapper(args): 
    custom_objects={"recall":recall,
                        "sensitivity":recall,
                        "specificity":specificity,
                        "fpr":fpr,
                        "fnr":fnr,
                        "precision":precision,
                        "f1":f1,
                        "ambig_binary_crossentropy":ambig_binary_crossentropy,
                        "ambig_mean_absolute_error":ambig_mean_absolute_error,
                        "ambig_mean_squared_error":ambig_mean_squared_error,
                        "MultichannelMultinomialMSE":MultichannelMultinomialMSE,
                        "RevCompConv1D":RevCompConv1D,
                        "BiasLayer":BiasLayer,
                        "MultichannelMultinomialNLL":MultichannelMultinomialNLL}
    get_custom_objects().update(custom_objects)
    if args.model_hdf5:
        model=load_model(args.model_hdf5)
        print(model.summary())
    else:
        model=model_from_json(open(args.json_string,"r").read())
        model=load_model_weights(args.weights,model)
    return model

def get_generator(args):
    gen=TiledbPredictGenerator(ref_fasta=args.ref_fasta,
                               batch_size=args.batch_size,
                               bed_regions_center=args.bed_regions_center,
                               bed_regions=args.bed_regions,
                               tdb_partition_thresh_for_upsample=None,
                               tdb_partition_attribute_for_upsample=None,
                               tdb_partition_datasets_for_upsample=None,
                               tdb_output_datasets=args.tdb_output_datasets,
                               tdb_input_datasets=args.tdb_input_datasets,
                               tdb_array=args.tdb_array,
                               chrom_sizes=args.chrom_sizes,
                               tdb_input_flank=args.tdb_input_flank,
                               tdb_input_source_attribute=args.tdb_input_source_attribute,
                               tdb_input_aggregation=args.tdb_input_aggregation,
                               tdb_input_transformation=args.tdb_input_transformation,
                               tdb_output_source_attribute=args.tdb_output_source_attribute,
                               tdb_output_flank=args.tdb_output_flank,
                               tdb_output_aggregation=args.tdb_output_aggregation,
                               tdb_output_transformation=args.tdb_output_transformation,
                               num_inputs=args.num_inputs,
                               num_outputs=args.num_outputs,
                               upsample_ratio=None,
                               tdb_ambig_attribute=None,
                               tdb_config=get_default_config(),
                               tdb_ctx=tiledb.Ctx(config=get_default_config()),
                               num_threads=args.num_threads)
    return gen

def get_interpretations(gen, model, count_explainer, prof_explainer,task_index,chipseq, tfchipseq):
    label_prof_dict={} 
    label_count_dict={} 
    pred_prof_dict={} 
    pred_count_dict={} 
    profile_shap_dict={}
    count_shap_dict={}
    seq_dict={}
    length_gen=len(gen)
    for i in range(length_gen): 
        print(str(i)+'/'+str(length_gen))
        X,y,coords=gen[i]
        coords=[[entry.decode('utf8')  for entry in coord] for coord in coords]
        preds=model.predict(X)

        pred_prob=get_probability_track_from_bpnet(preds[0][:,:,task_index])
        label_prob=get_probability_label_track(y[0][:,:,task_index])
        
        label_sum=y[1][:,task_index]
        pred_sum=preds[1][:,task_index]
        #chipseq will have control profile, control counts; ATAC/histone will not
        seq_input=X[0]
        if chipseq is True:
            control_profile=np.zeros_like(X[1])
            control_counts=np.zeros_like(X[2])
            count_explanations=count_explainer.shap_values([seq_input,control_counts])[0]
            profile_explanations=prof_explainer(seq_input, control_profile)
            print(profile_explanations.shape)
            profile_explanations=profile_explanations[task_index]
            count_explanations=count_explanations[task_index]       
        elif tfchipseq is True:
            control_profile=np.zeros_like(X[2])
            control_counts=np.zeros_like(X[1])
            count_explanations=count_explainer.shap_values([seq_input,control_counts])[0]
            profile_explanations=prof_explainer(seq_input, control_profile)
            profile_explanations=profile_explanations[0]
            count_explanations=count_explanations[0]       
        else:
            control_profile=None
            count_explanations=count_explainer.shap_values(X)[0]
            profile_explanations=prof_explainer(seq_input, control_profile)


        #store outputs in dictionary 
        for coord_index in range(len(coords)): 
            cur_coord=coords[coord_index][0:2]
            cur_coord[1]=int(cur_coord[1])
            cur_coord=tuple(cur_coord)
            label_prof_dict[cur_coord]=label_prob[coord_index]
            label_count_dict[cur_coord]=label_sum[coord_index]
            pred_prof_dict[cur_coord]=pred_prob[coord_index]
            pred_count_dict[cur_coord]=pred_sum[coord_index]
            profile_shap_dict[cur_coord]=profile_explanations[coord_index,:]
            count_shap_dict[cur_coord]=count_explanations[coord_index,:]
            seq_dict[cur_coord]=X[0][coord_index]
    return label_prof_dict, label_count_dict,pred_prof_dict,pred_count_dict, profile_shap_dict, count_shap_dict, seq_dict

def main():
    args=parse_args()
    gen=get_generator(args)
    print("created generator")
    #load the model
    model=load_model_wrapper(args)
    print("loaded model")    
    if args.chipseq is True:
        create_background_counts=create_background_counts_chip_1
        model_wrapper_for_counts=([model.input[0],model.input[2]],model.outputs[1][:,args.task_index:args.task_index+1])
        prof_explainer = create_explainer(model,ischip=args.chipseq,task_index=args.task_index)

    elif args.tfchipseq is True:
        create_background_counts=create_background_counts_chip_1
        model_wrapper_for_counts=([model.input[0],model.input[1]],model.outputs[0][:,args.task_index:args.task_index+1])
        prof_explainer = create_explainer_tf(model,ischip=args.tfchipseq,task_index=args.task_index)

    else:
        create_background_counts=create_background_atac 
        #oprint(model.outputs[1])
        model_wrapper_for_counts=(model.input, model.outputs[1][:,args.task_index:args.task_index+1])    
        prof_explainer = create_explainer(model,ischip=args.chipseq,task_index=args.task_index)

    count_explainer=shap.DeepExplainer(model_wrapper_for_counts,data=create_background_counts,combine_mult_and_diffref=combine_mult_and_diffref_1d)
    print("got count explainer") 
    print("got profile explainer")
    label_prof_dict, label_count_dict,pred_prof_dict,pred_count_dict, profile_shap_dict, count_shap_dict, seq_dict=get_interpretations(gen,model, count_explainer,prof_explainer,args.task_index,args.chipseq,args.tfchipseq)
    print("finished with interpretations")
    #save the dictionaries to disk! 
    
    outputs={} 
    outputs['label_prof']=label_prof_dict
    outputs['label_sum']=label_count_dict
    outputs['pred_prof']=pred_prof_dict
    outputs['pred_sum']=pred_count_dict
    outputs['profile_shap']=profile_shap_dict 
    outputs['count_shap']=count_shap_dict 
    outputs['seq']=seq_dict 
    pickle.dump(outputs,open(args.out_pickle, "wb" ) )


if __name__=="__main__":
    main()
    
