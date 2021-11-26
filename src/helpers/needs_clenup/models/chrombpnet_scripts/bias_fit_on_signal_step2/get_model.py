import importlib
import imp

from kerasAC.s3_sync import *
from kerasAC.custom_losses import *
from kerasAC.metrics import *
#from profile_bpnet_dnase_with_bias import attribution_prior_model
import keras

def get_w1_w0_training(args,train_generator):
    w1=args.w1
    w0=args.w0
    if (args.w1_w0_file is not None) and (args.w1_w0_file.startswith('s3://')):
        w1_w0_file=download_s3_file(args.w1_w0_file)
    else:
        w1_w0_file=args.w1_w0_file
    if (args.weighted==True and (w1==None or w0==None) ):
        if args.w1_w0_file==None:
            w1=train_generator.w1
            w0=train_generator.w0        
            assert args.save_w1_w0 !=None
            with open(args.save_w1_w0, 'w') as weight_file:
                for i in range(len(w1)):
                    weight_file.write(str(w1[i])+'\t'+str(w0[i])+'\n')
        else:
            w1_w0=np.loadtxt(args.w1_w0_file)
            w1=list(w1_w0[:,0])
            w0=list(w1_w0[:,1]) 
        print("got weights!")
    return w1,w0

def get_w1_w0_prediction(args):
    w1=None
    w0=None
    if args.w1_w0_file is not None:
        if args.w1_w0_file.startswith('s3://'):
            w1_w0_file=download_s3_file(args.w1_w0_file)
        else:
            w1_w0_file=args.w1_w0_file
        w1_w0=np.loadtxt(w1_w0_file)
        w1=w1_w0[:,0]
        w0=w1_w0[:,1]
    if args.w1!=None:
        w1=args.w1
    if args.w0!=None:
        w0=args.w0 
    return w1,w0

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


def get_model(args):
    from tensorflow.keras.utils import get_custom_objects
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
                    "MultichannelMultinomialNLL":MultichannelMultinomialNLL,        
                    "BiasLayer":BiasLayer}        
    w1,w0=get_w1_w0_prediction(args)
    if type(w1) in [np.ndarray, list]: 
        loss_function=get_weighted_binary_crossentropy(w0,w1)
        custom_objects["weighted_binary_crossentropy"]=loss_function
    get_custom_objects().update(custom_objects)
    
    if args.yaml!=None:
        from tensorflow.keras.models import model_from_yaml
        #load the model architecture from yaml
        if args.yaml.startswith('s3://'):
            yaml_string=read_s3_file_contents(args.yaml)
        else: 
            yaml_string=open(args.yaml,'r').read()
        model=model_from_yaml(yaml_string)#,custom_objects=custom_objects) 
    elif args.json!=None:
        from tensorflow.keras.models import model_from_json
        #load the model architecture from json
        if args.json.startswith('s3://'):
            json_string=read_s3_file_contents(args.json)
        else: 
            json_string=open(args.json,'r').read()
        model=model_from_json(json_string)#,custom_objects=custom_objects)
    elif args.load_model_hdf5!=None: 
        #load from the hdf5
        from tensorflow.keras.models import load_model
        if args.load_model_hdf5.startswith("s3://"):
            model_hdf5=download_s3_file(args.load_model_hdf5)
        else:
            model_hdf5=args.load_model_hdf5
        custom_objects['attribution_prior_model'] = attribution_prior_model
        get_custom_objects().update(custom_objects)
        model=load_model(model_hdf5)
    else:
        #initialize model from user-supplied architecture
        try:
            if (args.architecture_from_file!=None):
                architecture_module=imp.load_source('',args.architecture_from_file)
            else:
                architecture_module=importlib.import_module('kerasAC.architectures.'+args.architecture_spec)
        except:
            raise Exception("could not import requested architecture, is it installed in kerasAC/kerasAC/architectures? Is the file with the requested architecture specified correctly?")
        model=architecture_module.getModelGivenModelOptionsAndWeightInits(args)
    if args.num_gpus >1:
        try:
            model=multi_gpu_model(model,gpus=args.num_gpus)
            print("Training/predicting on" +str(args.num_gpus)+" GPU's. Set args.multi_gpu = False to avoid this") 
        except:
            print("failed to instantiate multi-gpu model, defaulting to single-gpu model")
    model=load_model_weights(args.weights,model)
    print("prepared the model")
    model.summary()
    return model

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
