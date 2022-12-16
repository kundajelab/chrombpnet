import argparse
from argparse import RawTextHelpFormatter


desc = """==================================================================================================
		Bias factorized, base-resolution deep learning models of chromatin accessibility reveal 
		cis-regulatory sequence syntax, transcription factor footprints and regulatory variants
		=================================================================================================="""

def read_parser():

        parser = argparse.ArgumentParser(description=desc,formatter_class=RawTextHelpFormatter)
        subparsers = parser.add_subparsers(help="Must be eithier 'pipeline', 'train', 'bias', 'pred_bw', 'contribs_bw', 'motif_discovery' ,'footprints', or 'snp_score'.", required=True, dest='cmd')
        pipeline_parser = subparsers.add_parser("pipeline", help="End-to-end pipline with qc and evaluation with bias factorized ChromBPNet model")
        train_parser = subparsers.add_parser("train", help="Train bias factorized ChromBPNet model")
        
        bias_parser_full = subparsers.add_parser("bias", help="Tools to train and qc bias model")
        bias_parser_sub = bias_parser_full.add_subparsers(help="Must be eithier 'pipeline'  'train', 'qc'.", required=True, dest='cmd_bias')
        bias_parser = bias_parser_sub.add_parser("pipeline", help="End-to-end pipline with qc and evaluation with bias model")
        bias_parser_train = bias_parser_sub.add_parser("train", help="Train bias model")
        bias_parser_qc = bias_parser_sub.add_parser("qc", help="Do quality checks on the bias model")
        
        preds_parser = subparsers.add_parser("pred_bw", help="Make model predictions on given regions and output to bigwig file")
        contribs_parser = subparsers.add_parser("contribs_bw", help="Get contribution score bigwigs")
        motifs_parser = subparsers.add_parser("modisco_motifs", help="Summarize motifs with TFModisco")
        #custom_preds_parser = subparsers.add_parser("pred_custom", help="Make model predictions on custom sequences and output to .h5 file")
        #custom_contribs_parser = subparsers.add_parser("contribs_custom", help="Get contribution on custom sequences and output to .h5 file")
        footprints_parser = subparsers.add_parser("footprints", help="Get marginal footprinting for given model and given motifs")
        variants_parser = subparsers.add_parser("variant_score", help="Score SNPs with model ")

        def general_training_args(required_train, optional_train):
		
        	required_train.add_argument('-g','--genome', required=True, type=str, help="reference genome fasta file")
       		required_train.add_argument('-c', '--chrom-sizes', type=str, required=True, help="Chrom sizes file")
        	group = required_train.add_mutually_exclusive_group(required=True)
       		group.add_argument('-ibam', '--input-bam-file', type=str, help="Input BAM file")
        	group.add_argument('-ifrag', '--input-fragment-file', type=str, help="Input fragment file")
        	group.add_argument('-itag', '--input-tagalign-file', type=str, help="Input tagAlign file")
        	required_train.add_argument('-o', '--output-dir', type=str, required=True, help="Output dir (path/to/output/dir)")
        	required_train.add_argument('-d', '--data-type', required=True, type=str, choices=['ATAC', 'DNASE'], help="assay type")
        	required_train.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
        	required_train.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
        	required_train.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
 
        	optional_train.add_argument("-oth", "--outlier_threshold", type=float, default=0.9999, help="threshold to use to filter outlies")
        	optional_train.add_argument('-ps', '--plus-shift', type=int, default=None, help="Plus strand shift applied to reads. Estimated if not specified")
        	optional_train.add_argument('-ms', '--minus-shift', type=int, default=None, help="Minus strand shift applied to reads. Estimated if not specified")
        	optional_train.add_argument('--ATAC-ref-path', type=str, default=None, help="Path to ATAC reference motifs (ATAC.ref.motifs.txt used by default)")
        	optional_train.add_argument('--DNASE-ref-path', type=str, default=None, help="Path to DNASE reference motifs (DNASE.ref.motifs.txt used by default)")
        	optional_train.add_argument('--num-samples', type=int, default=10000, help="Number of reads to sample from BAM/fragment/tagAlign file for shift estimation")
        	optional_train.add_argument("-il", "--inputlen", type=int, default=2114, required=False, help="Sequence input length")
        	optional_train.add_argument("-ol", "--outputlen", type=int, default=1000, required=False, help="Prediction output length")
        	optional_train.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
      	 	optional_train.add_argument("-e", "--epochs", type=int, default=50, help="Maximum epochs to train")
        	optional_train.add_argument("-es", "--early-stop", type=int, default=5, help="Early stop limit, corresponds to 'patience' in callback")
        	optional_train.add_argument("-l", "--learning-rate", type=float, default=0.001, help="Learning rate for model training")
        	optional_train.add_argument("--trackables",nargs="*",default=['logcount_predictions_loss', 'loss', 'logits_profile_predictions_loss', 'val_logcount_predictions_loss', 'val_loss', 'val_logits_profile_predictions_loss'], help="list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,val_logcount_predictions_loss,val_loss,val_profile_predictions_loss")
        	optional_train.add_argument("-a","--architecture_from_file",type=str,required=False, default=None, help="Model to use for training")
        	return required_train, optional_train
		
        # train chrombpnet arguments
		
        train_parser._action_groups.pop()
        required_main_parser = train_parser.add_argument_group('required arguments')
        optional_main_parser = train_parser.add_argument_group('optional arguments')
        
        required_main_parser,optional_main_parser = general_training_args(required_main_parser, optional_main_parser)
        
        required_main_parser.add_argument("-b", "--bias_model_path", type=str, required=True, help="Path for a pretrained bias model")

        optional_main_parser.add_argument("-sr", "--negative_sampling_ratio", type=float, default=0.1, help="Ratio of negatives to positive samples per epoch")
        optional_main_parser.add_argument("-fil", "--filters", type=int, default=512, help="Number of filters to use in chrombpnet mode")
        optional_main_parser.add_argument("-dil", "--n_dilation_layers", type=int, default=8, help="Number of dilation layers to use in chrombpnet model")
        optional_main_parser.add_argument("-j", "--max_jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_main_parser.add_argument("-bs", "--batch_size", type=int, default=64, help="batch size to use for model training")
 
         # chrombpnet pipeline arguments
		
        pipeline_parser._action_groups.pop()
        required_pipeline_parser = pipeline_parser.add_argument_group('required arguments')
        optional_pipeline_parser = pipeline_parser.add_argument_group('optional arguments')
        
        required_pipeline_parser,optional_pipeline_parser = general_training_args(required_pipeline_parser, optional_pipeline_parser)
        
        required_pipeline_parser.add_argument("-b", "--bias_model_path", type=str, required=True, help="Path for a pretrained bias model")

        optional_pipeline_parser.add_argument("-sr", "--negative_sampling_ratio", type=float, default=0.1, help="Ratio of negatives to positive samples per epoch")
        optional_pipeline_parser.add_argument("-fil", "--filters", type=int, default=512, help="Number of filters to use in chrombpnet mode")
        optional_pipeline_parser.add_argument("-dil", "--n_dilation_layers", type=int, default=8, help="Number of dilation layers to use in chrombpnet model")
        optional_pipeline_parser.add_argument("-j", "--max_jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_pipeline_parser.add_argument("-bs", "--batch_size", type=int, default=64, help="batch size to use for model training")
 
        # bias model pipeline arguments

        bias_parser._action_groups.pop()
        required_bias_parser = bias_parser.add_argument_group('required arguments')
        optional_bias_parser = bias_parser.add_argument_group('optional arguments')
        required_bias_parser,optional_bias_parser = general_training_args(required_bias_parser, optional_bias_parser)
        
        required_bias_parser.add_argument("-b", "--bias_threshold_factor", type=float, required=True, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions)")

        optional_bias_parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in chrombpnet mode")
        optional_bias_parser.add_argument("-dil", "--n_dilation_layers", type=int, default=4, help="Number of dilation layers to use in chrombpnet model")
        optional_bias_parser.add_argument("-j", "--max_jitter", type=int, default=0, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_bias_parser.add_argument("-bs", "--batch_size", type=int, default=64, help="batch size to use for model training")

       # bias model pipeline arguments

        bias_parser_train._action_groups.pop()
        required_biast_parser = bias_parser_train.add_argument_group('required arguments')
        optional_biast_parser = bias_parser_train.add_argument_group('optional arguments')
        required_biast_parser,optional_bias_parser = general_training_args(required_biast_parser, optional_biast_parser)
        
        required_biast_parser.add_argument("-b", "--bias_threshold_factor", type=float, required=True, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions)")

        optional_bias_parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in chrombpnet mode")
        optional_bias_parser.add_argument("-dil", "--n_dilation_layers", type=int, default=4, help="Number of dilation layers to use in chrombpnet model")
        optional_bias_parser.add_argument("-j", "--max_jitter", type=int, default=0, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_bias_parser.add_argument("-bs", "--batch_size", type=int, default=64, help="batch size to use for model training")
 

        # Make prediction bigwigs
        
        preds_parser._action_groups.pop()
        required_preds = preds_parser.add_argument_group('required arguments')
        optional_preds = preds_parser.add_argument_group('optional arguments')

        required_preds.add_argument("-bm", "--bias-model", type=str, required=False, help="Path to bias model h5 (one of -bm, -cm, -cmb reqd)")
        required_preds.add_argument("-cm", "--chrombpnet-model", type=str, required=False, help="Path to chrombpnet model h5 (one of -bm, -cm, -cmb reqd)")
        required_preds.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=False, help="Path to chrombpnet no bias model h5 (one of -bm, -cm, -cmb reqd)")
        required_preds.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of regions for prediction")
        required_preds.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_preds.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
        required_preds.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        optional_preds.add_argument("-os", "--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
        optional_preds.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        optional_preds.add_argument("-t", "--tqdm", type=int,default=1, help="Use tqdm. If yes then you need to have it installed.")
        optional_preds.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
        
        # Make contribution score bigwigs
        
        contribs_parser._action_groups.pop()
        required_contribs = contribs_parser.add_argument_group('required arguments')
        optional_contribs = contribs_parser.add_argument_group('optional arguments')

        required_contribs.add_argument("-m", "--model_h5", type=str, required=False, help="Path model .h5 file")
        required_contribs.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of regions for contribution score predictions")
        required_contribs.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_contribs.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
        required_contribs.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        optional_contribs.add_argument("-pc", "--profile_or_counts", nargs="+", type=str, default=["counts", "profile"], choices=["counts", "profile"],
                        help="use either counts or profile or both for running shap") 
        optional_contribs.add_argument("-os", "--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
        optional_contribs.add_argument("-t", "--tqdm", type=int,default=1, help="Use tqdm. If yes then you need to have it installed.")
        optional_contribs.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    
   
        # Get marginal footprints
        
        footprints_parser._action_groups.pop()
        required_ftps = footprints_parser.add_argument_group('required arguments')
        optional_ftps = footprints_parser.add_argument_group('optional arguments')

        required_ftps.add_argument("-m", "--model_h5", type=str, required=False, help="Path model .h5 file")
        required_ftps.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of regions for contribution score predictions")
        required_ftps.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_ftps.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
        required_ftps.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        required_ftps.add_argument("-pwm_f", "--motifs_to_pwm", type=str, required=True, help="Path to a TSV file containing motifs in first column and motif string to use for footprinting in second column")    
        
        optional_ftps.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        optional_ftps.add_argument("--ylim",default=None,type=tuple, required=False,help="lower and upper y-limits for plotting the motif footprint, in the form of a tuple i.e. \
        (0,0.8). If this is set to None, ylim will be autodetermined.")
  
        # Do variant scoring
        
        variants_parser._action_groups.pop()
        required_ves =  variants_parser.add_argument_group('required arguments')
        optional_ves =  variants_parser.add_argument_group('optional arguments')
        required_ves.add_argument("-snps", "--snp_data", type=str, required=True, help="Path to a tsv output with the following information in columns - chr, position to insert allele (0-based), ref allele, alt allele")
        required_ves.add_argument("-m", "--model_h5", type=str, required=False, help="Path model .h5 file")
        required_ves.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_ves.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
   
        optional_ves.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        optional_ves.add_argument("-dm","--debug_mode_on", type=int, default=0, help="Use this mode to print the flanks of first five SNP insert locations")
        
        
        # Run TF-Modisco
        
        motifs_parser._action_groups.pop()
        required_tfm =  motifs_parser.add_argument_group('required arguments')
        optional_tfm =  motifs_parser.add_argument_group('optional arguments')
        
        required_tfm.add_argument("-i", "--h5py", type=str, help="A legacy h5py file containing the one-hot encoded sequences and shap scores.")
        required_tfm.add_argument("-n", "--max_seqlets", type=int, required=True, help="The maximum number of seqlets per metacluster.")
        required_tfm.add_argument("-op", "--output-prefix", type=str, required=True, help="The path to the output file.")

        optional_tfm.add_argument("-l", "--n_leiden", type=int, default=2, help="The number of Leiden clusterings to perform with different random seeds.")
        optional_tfm.add_argument("-w", "--window", type=int, default=500, help="The window surrounding the peak center that will be considered for motif discovery.")
        optional_tfm.add_argument("-v", "--verbose", action="store_true", default=False, help="Controls the amount of output from the code.")

        
        # Pull the arguments
        
        args = parser.parse_args()
        print(args.cmd)
        if "cmd_bias" in args:
        	print(args.cmd_bias)
        return args

