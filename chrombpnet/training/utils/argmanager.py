import argparse

def update_data_args(parser):
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=False, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, default="None", help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, default="None" ,help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - see splits.py to set folds")


def update_train_args(parser):
    parser.add_argument("-e", "--epochs", type=int, default=50, help="Maximum epochs to train")
    parser.add_argument("-es", "--early-stop", type=int, default=5, help="Early stop limit, corresponds to 'patience' in callback")
    parser.add_argument("-bs", "--batch_size", type=int, default=64)
    parser.add_argument("-l", "--learning-rate", type=float, default=0.001)
    parser.add_argument("-pf", "--params", type=str, required=True, default=None)
    parser.add_argument("--trackables",nargs="*",default=['loss','val_loss'], help="list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,val_logcount_predictions_loss,val_loss,val_profile_predictions_loss")
    parser.add_argument("--optimizer", type=str, default="adam", choices=["adam", "muon"], help="adam, or muon for the dilated conv kernels (Adam at --learning-rate for the rest)")
    parser.add_argument("--muon-lr", type=float, default=None, help="Learning rate of the Muon-optimized conv kernels (default 2e-3)")
    parser.add_argument("--ema", action="store_true", help="Validate and save an exponential moving average of the weights")
    parser.add_argument("--lr-schedule", type=str, default="constant", choices=["constant", "cosine"], help="constant, or warmup + cosine decay over --epochs")
    parser.add_argument("--precision", type=str, default="default", choices=["default", "highest", "bf16"], help="default (TF32 on recent GPUs), highest (full float32) or bf16 (mixed bfloat16)")
    parser.add_argument("--device", type=str, default="auto", choices=["auto", "gpu", "cpu"], help="gpu fails if JAX has no GPU")

def update_model_args(parser):
    parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
    parser.add_argument("-a","--architecture_from_file",type=str,required=True, default=None, help="Model to use for training")

def fetch_train_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    update_train_args(parser)
    update_model_args(parser)
    args = parser.parse_args()

    assert((args.peaks.lower() != "none") or (args.nonpeaks.lower() != "none")) #Both peaks and nonpeaks are empty" 

    return args

def fetch_predict_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to model hdf5")
    parser.add_argument("-bs", "--batch_size", type=int, default=512)
    parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    args = parser.parse_args()

    assert((args.peaks.lower() != "none") or (args.nonpeaks.lower() != "none")) #Both peaks and nonpeaks are empty" 
    assert(args.inputlen % 2 ==0)
    assert(args.outputlen % 2 ==0)

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
