import argparse

def update_data_args(parser):
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=False, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-tdb_a", "--tdb_array", type=str, required=False, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-tdb_d", "--tdb_dataset", type=str, required=False, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-gen", "--generator", type=str, required=True, help="Generator to use - tiledb or batchgen")
    parser.add_argument("-fl", "--fold", type=str, required=True, help="Fold information - see splits.py to set folds")
    parser.add_argument("--trackables",nargs="*",default=['loss','val_loss'], help="list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,val_logcount_predictions_loss,val_loss,val_profile_predictions_loss")


def update_train_args(parser):
    parser.add_argument("-e", "--epochs", type=int, default=50, help="Maximum epochs to train")
    parser.add_argument("-es", "--early-stop", type=int, default=5, help="Early stop limit, corresponds to 'patience' in callback")
    parser.add_argument("-bs", "--batch-size", type=int, default=64)
    parser.add_argument("-l", "--learning-rate", type=float, default=0.001)
    parser.add_argument("-pf", "--params", type=str, required=True, default=None)

def update_model_args(parser):
    parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("-a","--architecture_from_file",type=str,required=True, default=None, help="Model to use for training")

def fetch_train_chrombpnet_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    update_train_args(parser)
    update_model_args(parser)

    # additional arguments for chrombpnet
    parser.add_argument("-j", "--max-jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet")    
    parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=1.0, help="Ratio of negative to positive samples per epoch")
    parser.set_defaults(negative_sampling=True)

    args = parser.parse_args()
    return args

def fetch_predict_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to model hdf5")
    parser.add_argument("-bs", "--batch-size", type=int, default=512)
    parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    args = parser.parse_args()
    return args


def fetch_modisco_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--scores-prefix", type=str, required=True, help="Prefix to counts/profile h5 files. Will use prefix.{profile,counts}_scores.h5")
    parser.add_argument("-p","--profile-or-counts", type=str, required=True, help="Scoring method to use, profile or counts scores")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("-c", "--crop", type=int, default=500, help="Crop scores to this width from the center for each example")
    parser.add_argument("-m", "--max-seqlets", type=int, default=50000, help="Max number of seqlets per metacluster for modisco")

    args = parser.parse_args()
    return args


def update_tiledb_args():
    # freeze tiledb to use only one dataset
    parser=argparse.ArgumentParser()
    
    parser.add_argument("--num_inputs",type=int, default=1)
    parser.add_argument("--num_outputs",type=int, default=2)

    tiledbgroup=parser.add_argument_group('tiledb')
    tiledbgroup.add_argument("-tdb", "--tiledb_array", type=str, required=False, help="Needed if tiledb generator is used")
    tiledbgroup.add_argument("--tdb_output_datasets",nargs="+", required=False, default=None, help="dataset column from db_ingest; comma separated across channels; space separated across outputs")
    tiledbgroup.add_argument("--tdb_output_source_attribute", required=False, nargs="+",help="tiledb attribute for use in label generation i.e. fc_bigwig")        

    tiledbgroup.add_argument("--tdb_input_datasets",nargs="+",default=None,help="dataset column from db_ingest; comma separated across channels; space separated across inputs")
    tiledbgroup.add_argument("--tdb_input_source_attribute",nargs="+",help="attribute to use for generating model input, or 'seq' for one-hot-encoded sequence")

    arch_params.add_argument("--model_params",type=str,default=None,help="2-column file with param name in column 1 and param value in column 2")
    arch_params.add_argument("--datasets",nargs="*",default=None,help="list of tasks to train on, by name")

    batch_params=parser.add_argument_group("batch_params")
    batch_params.add_argument("--batch_size",type=int,default=1000)
    batch_params.add_argument("--revcomp",action="store_true")

    parallelization_params=parser.add_argument_group("parallelization")
    parallelization_params.add_argument("--threads",type=int,default=1)
    parallelization_params.add_argument("--max_queue_size",type=int,default=100)

    vis_params.add_argument("--trackables",nargs="*",default=['loss','val_loss'], help="list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,val_logcount_predictions_loss,val_loss,val_profile_predictions_loss")
    return parser.parse_args()
