from __future__ import division, print_function, absolute_import
import importlib.util
import re
import sys
import jax
import keras
import chrombpnet.training.utils.argmanager as argmanager
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.callbacks as callbacks
import chrombpnet.training.data_generators.initializers as initializers
from chrombpnet.training import runtime
from chrombpnet.training.utils.model_io import require_jax_backend
import pandas as pd
import os
import json
import numpy as np

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
os.environ['PYTHONHASHSEED'] = '0'

def get_model(args, parameters):
    """
    Read a model definition from a python file. This function can be used to read any model architecture that takes sequence as input
    and outputs a two task model.  Task one to predict the probability distribution of a profile and task two to predict the total counts in a profile.
    Look at .py models in src/training/models/ for examples. I will try to provide a dummy model as example here - for later.
    The files should have the following two functions - getModelGivenModelOptionsAndWeightInits and save_model_without_bias
    """
    # load under a real module name (not ''), so objects defined in the file get a usable __module__
    stem = re.sub(r"\W", "_", os.path.splitext(os.path.basename(args.architecture_from_file))[0])
    spec = importlib.util.spec_from_file_location("chrombpnet_architecture_" + stem, args.architecture_from_file)
    architecture_module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = architecture_module
    spec.loader.exec_module(architecture_module)
    model=architecture_module.getModelGivenModelOptionsAndWeightInits(args, parameters)
    print("got the model")
    return model, architecture_module

class Float32ModelCheckpoint(keras.callbacks.ModelCheckpoint):
    """ModelCheckpoint that writes the model with float32 policies (used under --precision bf16), so that the
    best-so-far model left by an interrupted run loads and runs in float32 like the final one. The layers go back
    to their training policies right after each save."""

    def _save_model(self, epoch, batch, logs):
        # every ModelCheckpoint save goes through _save_model (tests/test_training.py checks the interrupted case)
        with runtime.float32_policy(self.model):
            super()._save_model(epoch, batch, logs)

def fit_and_evaluate(model,train_gen,valid_gen,args,architecture_module):
    model_output_path_h5_name=args.output_prefix+".h5"
    model_output_path_logs_name=args.output_prefix+".log"

    checkpoint_class = Float32ModelCheckpoint if runtime.bf16_active() else keras.callbacks.ModelCheckpoint
    checkpointer = checkpoint_class(filepath=model_output_path_h5_name, monitor="val_loss", mode="min",  verbose=1, save_best_only=True)
    # Keras 3 restores the best weights at the end of training even when --epochs is reached without an early stop
    earlystopper = keras.callbacks.EarlyStopping(monitor='val_loss', mode="min", patience=args.early_stop, verbose=1, restore_best_weights=True)
    history= callbacks.LossHistory(model_output_path_logs_name+".batch",args.trackables)
    csvlogger = keras.callbacks.CSVLogger(model_output_path_logs_name, append=False)
    #reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.4, patience=args.early_stop-2, min_lr=0.00000001)
    cur_callbacks=[checkpointer,earlystopper,csvlogger,history]
    if getattr(args, "ema", False):
        # validate on the EMA weights and put them in the model at every epoch end, before the checkpoint and the
        # early stopping callbacks see the model, so the saved models hold EMA weights
        cur_callbacks.insert(0, keras.callbacks.SwapEMAWeights(swap_on_epoch=True))

    model.fit(train_gen,
              validation_data=valid_gen,
              epochs=args.epochs,
              verbose=1,
              callbacks=cur_callbacks)

    if runtime.bf16_active():
        # trained with mixed_bfloat16: save float32 models (the weights are float32 already) so that predict,
        # footprints and interpret run them in full precision
        runtime.set_float32_policy(model)

    print('save model') 
    model.save(model_output_path_h5_name)

    architecture_module.save_model_without_bias(model, args.output_prefix)


def get_model_param_dict(args):
    '''
    param_file is a TSV file with 2 columns -- param name in column 1, and param value in column 2
    You can pass model specfic parameters to design your own model with this.
    '''
    params={}
    for line in open(args.params,'r').read().strip().split('\n'):
        tokens=line.split('\t')
        params[tokens[0]]=tokens[1]

    assert("counts_loss_weight" in params.keys()) # missing counts loss weight to use
    assert("filters" in params.keys()) # filters to use for the model not provided
    assert("n_dil_layers" in params.keys()) # n_dil_layers to use for the model not provided
    assert("inputlen" in params.keys()) # inputlen to use for the model not provided
    assert("outputlen" in params.keys()) # outputlen to use for the model not provided
    assert("negative_sampling_ratio" in params.keys()) # negative_sampling_ratio to use for the model not provided
    assert("max_jitter" in params.keys()) # max_jitter to use for the model not provided
    assert(args.chr_fold_path==params["chr_fold_path"]) # the parameters were generated on a different folds compared to the given fold

    assert(int(params["inputlen"])%2==0)
    assert(int(params["outputlen"])%2==0)

    return params 

def json_safe_args(args):
    """args as a JSON-serializable dict: primitives (and lists of them) as they are, anything else as str."""
    primitive = (str, int, float, bool, type(None))
    out = {}
    for key, value in vars(args).items():
        if isinstance(value, (list, tuple)) and all(isinstance(v, primitive) for v in value):
            out[key] = list(value)
        elif isinstance(value, primitive):
            out[key] = value
        else:
            out[key] = str(value)
    return out

def main(args):

    require_jax_backend()
    # training needs the backends now: create and log them even for --device auto
    runtime.assert_gpu_if_requested(getattr(args, "device", None) or "auto", initialize=True)
    # --precision bf16 / highest set a process-wide dtype policy / matmul precision; restore them afterwards so
    # later pipeline steps (predict, interpret) build, load and run models as before
    previous_policy = keras.config.dtype_policy()
    previous_matmul_precision = jax.config.jax_default_matmul_precision
    runtime.configure_precision(getattr(args, "precision", None) or "default")
    try:
        # read tab-seperated parameters file
        parameters = get_model_param_dict(args)
        print(parameters)
        # seeds python random (initial weights, batch order), numpy (jitter, revcomp, negatives) and keras
        keras.utils.set_random_seed(args.seed)

        train_generator = None
        if getattr(args, "lr_schedule", None) == "cosine" and not getattr(args, "total_steps", None):
            # the schedule is part of the optimizer, which the architecture file builds: it needs the step budget
            train_generator = initializers.initialize_generators(args, "train", parameters, return_coords=False)
            args.total_steps = len(train_generator) * args.epochs
            numpy_state = np.random.get_state()

        # get model architecture to load
        model, architecture_module=get_model(args, parameters)

        # initialize generators to load data
        if train_generator is None:
            train_generator = initializers.initialize_generators(args, "train", parameters, return_coords=False)
        else:
            # the architecture file reseeds numpy: continue the stream the training generator started, so the
            # next epochs do not repeat its jitter / revcomp draws
            np.random.set_state(numpy_state)
        valid_generator = initializers.initialize_generators(args, "valid", parameters, return_coords=False)

        # train the model using the generators
        fit_and_evaluate(model, train_generator, valid_generator, args, architecture_module)
        run_info = runtime.runtime_info()
    finally:
        keras.config.set_dtype_policy(previous_policy)
        jax.config.update("jax_default_matmul_precision", previous_matmul_precision)

    # store arguments and and parameters to checkpoint
    with open(args.output_prefix+'.args.json', 'w') as fp:
        json.dump(dict(json_safe_args(args), runtime=run_info), fp,  indent=4)
    #with open(args.output_prefix+'.params.json', 'w') as fp:
    #    json.dump(parameters, fp,  indent=4)


if __name__=="__main__":
    # read arguments
    args=argmanager.fetch_train_args()
    main(args)

