import argparse

def update_data_args(parser):
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-o", "--output-prefix", type=str, required=True, help="Output prefix")

def update_train_args(parser):
    parser.add_argument("-tc", "--test-chr", nargs="+", default=["chr1"], help="Test chromosome/s")
    parser.add_argument("-vc", "--val-chr", nargs="+", default=["chr8", "chr10"], help="Validattion chromosome/s")
    parser.add_argument("-w", "--counts-weight", type=float, default=0.1, help="Counts MSE loss is assigned a weight of mean(total_counts_per_peak) * counts_weight relative to multinomial loss")
    parser.add_argument("-e", "--epochs", type=int, default=50, help="Maximum epochs to train")
    parser.add_argument("-es", "--early-stop", type=int, default=5, help="Early stop limit, corresponds to 'patience' in callback")
    parser.add_argument("-bs", "--batch-size", type=int, default=64)
    parser.add_argument("-l", "--learning-rate", type=float, default=0.001)

def update_model_args(parser):
    parser.add_argument("-il", "--inputlen", type=int, default=2114, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, default=1000, help="Prediction output length")
    parser.add_argument("-f", "--filters", type=int, default=512, help="Number of filters for BPNet")
    parser.add_argument("-nd", "--ndil", type=int, default=8, help="Number of dilated convs for BPNet")

def fetch_train_bias_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    update_train_args(parser)
    update_model_args(parser)
    parser.add_argument("-t", "--threshold-factor", type=float, default=0.5, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions)")
    parser.add_argument("-j", "--max-jitter", type=int, default=100, help="Maximum jitter applied on either side of region (default 100 for bias model)")
    args = parser.parse_args()
    return args

def fetch_train_chrombpnet_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    update_train_args(parser)
    update_model_args(parser)

    # additional arguments for chrombpnet
    parser.add_argument("-j", "--max-jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet")
    parser.add_argument("-bm", "--bias-model", type=str, required=True, help="Path to bias model hdf5")
    parser.add_argument("--no-adjust-bias-cts", dest='adjust_bias_cts', action='store_false', 
                        help="""NOT RECOMMENDED. By default the code adjusts bias predicted counts to 
                              account for different sequencing depths between bias model source bigwig
                              and current bigwig. You may choose to not do this if the bias model was 
                              trained on the same data as the current data, but even then it is recommended 
                              to adjust counts.""")
    parser.set_defaults(adjust_bias_cts=True)
    
    parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=0.1, help="Ratio of negative to positive samples per epoch")
    parser.add_argument("--no-negative-sampling", dest='negative_sampling', action='store_false', help="Include this flag if sampling of negatives is not desired. This will use all negatives supplied every epoch.")
    parser.set_defaults(negative_sampling=True)

    args = parser.parse_args()
    return args

def fetch_metrics_args():
    parser = argparse.ArgumentParser()
    update_data_args(parser)
    parser.add_argument("-tc", "--test-chr", nargs="+", required=True, help="Test chromosome/s, e.g. -tc chr1, space separated for more than 1")
    parser.add_argument("-bm", "--bias-model", type=str, required=True, help="Path to bias model hdf5")
    parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model hdf5")
    parser.add_argument("-bs", "--batch-size", type=int, default=512)

    args = parser.parse_args()
    return args

def fetch_interpret_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-m", "--model", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-o", "--output-prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    parser.add_argument("-p", "--profile-or-counts", nargs="+", type=str, default=["counts", "profile"], choices=["counts", "profile"],
                        help="use either counts or profile or both for running shap")

    args = parser.parse_args()
    return args

def fetch_modisco_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--scores-prefix", type=str, required=True, help="Prefix to counts/profile h5 files. Will use prefix.{profile,counts}_scores.h5")
    parser.add_argument("-p","--profile-or-counts", type=str, required=True, help="Scoring method to use, profile or counts scores")
    parser.add_argument("-o", "--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("-c", "--crop", type=int, default=500, help="Crop scores to this width from the center for each example")
    parser.add_argument("-m", "--max-seqlets", type=int, default=1000000, help="Max number of seqlets per metacluster for modisco")

    args = parser.parse_args()
    return args
