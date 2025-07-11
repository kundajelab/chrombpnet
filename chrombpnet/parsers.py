import argparse
from argparse import RawTextHelpFormatter


desc = """==================================================================================================
		Bias factorized, base-resolution deep learning models of chromatin accessibility reveal 
		cis-regulatory sequence syntax, transcription factor footprints and regulatory variants
		=================================================================================================="""

def read_parser():

        parser = argparse.ArgumentParser(description=desc,formatter_class=RawTextHelpFormatter)
        subparsers = parser.add_subparsers(help="Must be eithier 'pipeline', 'train', 'qc', 'bias', 'prep', 'pred_bw', 'contribs_bw', 'modisco_motifs' ,'footprints', or 'snp_score'.", required=True, dest='cmd')
        
        # main parsers
        
        pipeline_parser = subparsers.add_parser("pipeline", help="End-to-end pipline with train, quality check and test for bias factorized ChromBPNet model")
        train_parser = subparsers.add_parser("train", help="Train bias factorized ChromBPNet model")
        qc_parser = subparsers.add_parser("qc", help="Do quality checks and get test metrics for bias factorized ChromBPNet model")

        # bias parsers
        
        bias_parser_full = subparsers.add_parser("bias", help="Tools to train, quality check and test bias model")
        bias_parser_sub = bias_parser_full.add_subparsers(help="Must be eithier 'pipeline'  'train', 'qc'.", required=True, dest='cmd_bias')
        bias_parser = bias_parser_sub.add_parser("pipeline", help="End-to-end pipline with train, quality check and test for bias model")
        bias_parser_train = bias_parser_sub.add_parser("train", help="Train bias model")
        bias_parser_qc = bias_parser_sub.add_parser("qc", help="Do quality checks and get test metrics for the bias model")
        
        # helper parsers

        prep_parser_full = subparsers.add_parser("prep", help="Tools to generate preprocessing data for chrombpnet")
        prep_parser_sub = prep_parser_full.add_subparsers(help="Must be eithier 'nonpeaks' or 'splits'.", required=True, dest='cmd_prep')
        nonpeaks_parser = prep_parser_sub.add_parser("nonpeaks", help="Generate non-peak background regions given peaks")
        splits_parser = prep_parser_sub.add_parser("splits", help="Generate chromosome splits")   
        
        # downstream tool parsers
        preds_parser = subparsers.add_parser("pred_bw", help="Get model prediction bigwigs (Metrics calculated if observed bigwig provided)")
        contribs_parser = subparsers.add_parser("contribs_bw", help="Get contribution score bigwigs")
        #motifs_parser = subparsers.add_parser("modisco_motifs", help="(Will soon be deprecated: use modisco motifs from tfmodisco lite) Summarize motifs from contribution scores with TFModisco")
        #custom_preds_parser = subparsers.add_parser("pred_custom", help="Make model predictions on custom sequences and output to .h5 file")
        #custom_contribs_parser = subparsers.add_parser("contribs_custom", help="Get contribution on custom sequences and output to .h5 file")
        footprints_parser = subparsers.add_parser("footprints", help="Get marginal footprinting for given model and given motifs")
        #variants_parser = subparsers.add_parser("snp_score", help="Score SNPs with model")

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
        	required_train.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
 
        	optional_train.add_argument("-oth", "--outlier-threshold", type=float, default=0.9999, help="threshold to use to filter outlies")
        	#optional_train.add_argument('-ps', '--plus-shift', type=int, default=None, help="Plus strand shift applied to reads. Estimated if not specified")
        	#optional_train.add_argument('-ms', '--minus-shift', type=int, default=None, help="Minus strand shift applied to reads. Estimated if not specified")
        	optional_train.add_argument('--ATAC-ref-path', type=str, default=None, help="Path to ATAC reference motifs (ATAC.ref.motifs.txt used by default)")
        	optional_train.add_argument('--DNASE-ref-path', type=str, default=None, help="Path to DNASE reference motifs (DNASE.ref.motifs.txt used by default)")
        	optional_train.add_argument('--num-samples', type=int, default=10000, help="Number of reads to sample from BAM/fragment/tagAlign file for shift estimation")
        	optional_train.add_argument("-il", "--inputlen", type=int, default=2114, required=False, help="Sequence input length")
        	optional_train.add_argument("-ol", "--outputlen", type=int, default=1000, required=False, help="Prediction output length")
        	optional_train.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for model training")
      	 	optional_train.add_argument("-e", "--epochs", type=int, default=50, help="Maximum epochs to train")
        	optional_train.add_argument("-es", "--early-stop", type=int, default=5, help="Early stop limit, corresponds to 'patience' in callback")
        	optional_train.add_argument("-l", "--learning-rate", type=float, default=0.001, help="Learning rate for model training")
        	optional_train.add_argument("-track","--trackables",nargs="*",default=['logcount_predictions_loss', 'loss', 'logits_profile_predictions_loss', 'val_logcount_predictions_loss', 'val_loss', 'val_logits_profile_predictions_loss'], help="list of things to track per batch, such as logcount_predictions_loss,loss,profile_predictions_loss,val_logcount_predictions_loss,val_loss,val_profile_predictions_loss")
        	optional_train.add_argument("-a","--architecture-from-file",type=str,required=False, default=None, help="Model to use for training")
        	optional_train.add_argument("-fp","--file-prefix",type=str,required=False, default=None, help="File prefix for output to use. All the files will be prefixed with this string if provided.")
        	optional_train.add_argument('-hp', '--html-prefix', required=False, default="./", help="The html prefix to use for the html file output.")
        	optional_train.add_argument('--bsort', required=False, default=False, action='store_true', help="In prpeprocess, by deafult we sort bam using unix sort but sometimes LC collate can cause issues, so this can be set to use betools sort which works well but is memory intensive..")
        	optional_train.add_argument('--tmpdir', required=False, default=None, type=str, help="temp dir for unix sort")
        	optional_train.add_argument('--no-st', required=False, default=False, action='store_true', help="Dont do streaming  and filtering in preprocessing (short chromosome contrigs not in reference fasta are not removed)")

        	return required_train, optional_train

        # Generate non-peak regions from peak-regions
        
        nonpeaks_parser._action_groups.pop()
        required_nonpeaks_parser = nonpeaks_parser.add_argument_group('required arguments')
        optional_nonpeaks_parser = nonpeaks_parser.add_argument_group('optional arguments')
        
        required_nonpeaks_parser.add_argument("-g","--genome", required=True, help="reference genome file")
        required_nonpeaks_parser.add_argument("-o","--output-prefix", required=True, help="output BED file prefix to store the gc content of binned genome. suffix .bed will be appended by the code. If the prefix contains a directory path make sure it exists.")
        required_nonpeaks_parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
        required_nonpeaks_parser.add_argument('-c', '--chrom-sizes', type=str, required=True, help="Chrom sizes file")
        required_nonpeaks_parser.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")

        optional_nonpeaks_parser.add_argument("-il","--inputlen", type=int,default=2114, help="inputlen to use to make bins and find gc content")
        optional_nonpeaks_parser.add_argument("-st","--stride", type=int,default=1000, help="stride to use for shifting the bins")
        optional_nonpeaks_parser.add_argument("-npr", "--neg-to-pos-ratio-train", type=int, default=2, help="Ratio of negatives to positives to sample in training set (test set always has 1:1 positive to negatives ratio")
        optional_nonpeaks_parser.add_argument("-br", "--blacklist-regions", type=str, required=False, default=None, help="TSV file with 3 columns - chr, start, end")
        optional_nonpeaks_parser.add_argument("-s", "--seed", type=int, default=1234, help="seed to use for generating nonpeaks")

        # Generate splits
        splits_parser._action_groups.pop()
        required_splits_parser = splits_parser.add_argument_group('required arguments')
        
        required_splits_parser.add_argument("-op", "--output_prefix", type=str, required=True, help="Path prefix to store the fold information (appended with .json)")
        required_splits_parser.add_argument("-c", "--chrom-sizes", type=str, required=True, help="TSV file with chromosome sizes. All chromosomes from the first column of chrom sizes file are used")
        required_splits_parser.add_argument("-tcr", "--test-chroms", nargs="*", type=str, required=True, help="Chromosomes to use for test")
        required_splits_parser.add_argument("-vcr", "--valid-chroms", nargs="*", type=str, required=True, help="Chromosomes to use for validation")

        # train chrombpnet arguments
		
        train_parser._action_groups.pop()
        required_main_parser = train_parser.add_argument_group('required arguments')
        optional_main_parser = train_parser.add_argument_group('optional arguments')
        
        required_main_parser,optional_main_parser = general_training_args(required_main_parser, optional_main_parser)
        
        required_main_parser.add_argument("-b", "--bias-model-path", type=str, required=True, help="Path for a pretrained bias model")

        optional_main_parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=0.1, help="Ratio of negatives to positive samples per epoch")
        optional_main_parser.add_argument("-fil", "--filters", type=int, default=512, help="Number of filters to use in chrombpnet mode")
        optional_main_parser.add_argument("-dil", "--n-dilation-layers", type=int, default=8, help="Number of dilation layers to use in chrombpnet model")
        optional_main_parser.add_argument("-j", "--max-jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_main_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")
 
         # chrombpnet pipeline arguments
		
        pipeline_parser._action_groups.pop()
        required_pipeline_parser = pipeline_parser.add_argument_group('required arguments')
        optional_pipeline_parser = pipeline_parser.add_argument_group('optional arguments')
        
        required_pipeline_parser,optional_pipeline_parser = general_training_args(required_pipeline_parser, optional_pipeline_parser)
        
        required_pipeline_parser.add_argument("-b", "--bias-model-path", type=str, required=True, help="Path for a pretrained bias model")

        optional_pipeline_parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=0.1, help="Ratio of negatives to positive samples per epoch")
        optional_pipeline_parser.add_argument("-fil", "--filters", type=int, default=512, help="Number of filters to use in chrombpnet mode")
        optional_pipeline_parser.add_argument("-dil", "--n-dilation-layers", type=int, default=8, help="Number of dilation layers to use in chrombpnet model")
        optional_pipeline_parser.add_argument("-j", "--max-jitter", type=int, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_pipeline_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")
 
        # chrombpnet model qc arguments

        qc_parser._action_groups.pop()
        required_qc_parser = qc_parser.add_argument_group('required arguments')
        optional_qc_parser = qc_parser.add_argument_group('optional arguments')

        required_qc_parser.add_argument("-bw", "--bigwig", type=str, required=True, help="Input bigwig file of observed data")
        required_qc_parser.add_argument("-cm", "--chrombpnet-model", type=str, required=True, help="Path to chrombpnet model h5")
        required_qc_parser.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=True, help="Path to chrombpnet nobias model h5")
        required_qc_parser.add_argument('-g','--genome', required=True, type=str, help="reference genome fasta file")
        required_qc_parser.add_argument('-c', '--chrom-sizes', type=str, required=True, help="Chrom sizes file")
        required_qc_parser.add_argument('-o', '--output-dir', type=str, required=True, help="Output dir (path/to/output/dir)")
        required_qc_parser.add_argument('-d', '--data-type', required=True, type=str, choices=['ATAC', 'DNASE'], help="assay type")
        required_qc_parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
        required_qc_parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
        required_qc_parser.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes") 
        optional_qc_parser.add_argument("-fp","--file-prefix",type=str,required=False, default=None, help="File prefix for output to use. All the files will be prefixed with this string if provided.")
        optional_qc_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")
        optional_qc_parser.add_argument('-hp', '--html-prefix', required=False, default="./", help="The html prefix to use for the html file output.")
 

        # bias model pipeline arguments

        bias_parser._action_groups.pop()
        required_bias_parser = bias_parser.add_argument_group('required arguments')
        optional_bias_parser = bias_parser.add_argument_group('optional arguments')
        required_bias_parser,optional_bias_parser = general_training_args(required_bias_parser, optional_bias_parser)
        
        required_bias_parser.add_argument("-b", "--bias-threshold-factor", type=float, required=True, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions). Recommended start value 0.5 for ATAC and 0.8 for DNase.")

        optional_bias_parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in chrombpnet mode")
        optional_bias_parser.add_argument("-dil", "--n-dilation-layers", type=int, default=4, help="Number of dilation layers to use in chrombpnet model")
        optional_bias_parser.add_argument("-j", "--max-jitter", type=int, default=0, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_bias_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")

       # bias model training arguments

        bias_parser_train._action_groups.pop()
        required_biast_parser = bias_parser_train.add_argument_group('required arguments')
        optional_biast_parser = bias_parser_train.add_argument_group('optional arguments')
        required_biast_parser,optional_biast_parser = general_training_args(required_biast_parser, optional_biast_parser)
        
        required_biast_parser.add_argument("-b", "--bias-threshold-factor", type=float, required=True, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions). Recommended start value 0.5 for ATAC and 0.8 for DNas")

        optional_biast_parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in chrombpnet mode")
        optional_biast_parser.add_argument("-dil", "--n-dilation-layers", type=int, default=4, help="Number of dilation layers to use in chrombpnet model")
        optional_biast_parser.add_argument("-j", "--max-jitter", type=int, default=0, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
        optional_biast_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")
 
       # bias model qc arguments

        bias_parser_qc._action_groups.pop()
        required_bqc_parser = bias_parser_qc.add_argument_group('required arguments')
        optional_bqc_parser = bias_parser_qc.add_argument_group('optional arguments')

        required_bqc_parser.add_argument("-bw", "--bigwig", type=str, required=True, help="Input bigwig file of observed data")
        required_bqc_parser.add_argument("-bm", "--bias-model", type=str, required=True, help="Path to bias model .h5 file")
        required_bqc_parser.add_argument('-g','--genome', required=True, type=str, help="reference genome fasta file")
        required_bqc_parser.add_argument('-c', '--chrom-sizes', type=str, required=True, help="Chrom sizes file")
        required_bqc_parser.add_argument('-o', '--output-dir', type=str, required=True, help="Output dir (path/to/output/dir)")
        required_bqc_parser.add_argument('-d', '--data-type', required=True, type=str, choices=['ATAC', 'DNASE'], help="assay type")
        required_bqc_parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
        required_bqc_parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
        required_bqc_parser.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes") 
        optional_bqc_parser.add_argument("-fp","--file-prefix",type=str,required=False, default=None, help="File prefix for output to use. All the files will be prefixed with this string if provided.")
        optional_bqc_parser.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for model training")
        optional_bqc_parser.add_argument('-hp', '--html-prefix', required=False, default="./", help="The html prefix to use for the html file output.")
 

        # Make prediction bigwigs
        
        preds_parser._action_groups.pop()
        required_preds = preds_parser.add_argument_group('required arguments')
        optional_preds = preds_parser.add_argument_group('optional arguments')

        required_preds.add_argument("-bm", "--bias-model", type=str, required=False, help="Path to bias model h5 (atleast one of -bm, -cm, -cmb  is reqd)")
        required_preds.add_argument("-cm", "--chrombpnet-model", type=str, required=False, help="Path to chrombpnet model h5 (atleast one of -bm, -cm, -cmb is reqd)")
        required_preds.add_argument("-cmb", "--chrombpnet-model-nb", type=str, required=False, help="Path to chrombpnet no bias model h5 (atleast one of -bm, -cm, -cmb  is reqd)")
        required_preds.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of regions for prediction")
        required_preds.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_preds.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
        required_preds.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        optional_preds.add_argument("-os", "--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
        optional_preds.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        optional_preds.add_argument("-t", "--tqdm", type=int,default=1, help="Use tqdm. If yes then you need to have it installed.")
        optional_preds.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
        optional_preds.add_argument("-bw", "--bigwig", type=str, default=None, help="If provided .h5 with predictions are output along with calculated metrics considering bigwig as groundtruth.")
       
        # Make contribution score bigwigs
        
        contribs_parser._action_groups.pop()
        required_contribs = contribs_parser.add_argument_group('required arguments')
        optional_contribs = contribs_parser.add_argument_group('optional arguments')

        required_contribs.add_argument("-m", "--model-h5", type=str, required=True, help="Path model .h5 file")
        required_contribs.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of regions for contribution score predictions")
        required_contribs.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_contribs.add_argument("-c", "--chrom-sizes", type=str, required=True, help="Chromosome sizes 2 column tab-separated file")
        required_contribs.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        optional_contribs.add_argument("-pc", "--profile-or-counts", nargs="+", type=str, default=["counts", "profile"], choices=["counts", "profile"],
                        help="use either counts or profile or both for running shap") 
        optional_contribs.add_argument("-os", "--output-prefix-stats", type=str, default=None, required=False, help="Output stats on bigwig")
        optional_contribs.add_argument("-t", "--tqdm", type=int,default=1, help="Use tqdm. If yes then you need to have it installed.")
        optional_contribs.add_argument("-d", "--debug-chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    
   
        # Get marginal footprints
        
        footprints_parser._action_groups.pop()
        required_ftps = footprints_parser.add_argument_group('required arguments')
        optional_ftps = footprints_parser.add_argument_group('optional arguments')

        required_ftps.add_argument("-m", "--model-h5", type=str, required=True, help="Path model .h5 file")
        required_ftps.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of non-peak regions")
        required_ftps.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        required_ftps.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
        required_ftps.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
        required_ftps.add_argument("-pwm_f", "--motifs-to-pwm", type=str, required=True, help="Path to a TSV file containing motifs in first column and motif string to use for footprinting in second column")    
        
        optional_ftps.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        optional_ftps.add_argument("--ylim",default=None,type=tuple, required=False,help="lower and upper y-limits for plotting the motif footprint, in the form of a tuple i.e. \
        (0,0.8). If this is set to None, ylim will be autodetermined.")
  
        # Do variant scoring
        
        #variants_parser._action_groups.pop()
        #required_ves =  variants_parser.add_argument_group('required arguments')
        #optional_ves =  variants_parser.add_argument_group('optional arguments')
        #required_ves.add_argument("-snps", "--snp-data", type=str, required=True, help="Path to a tsv output with the following information in columns - chr, position to insert allele (0-based), ref allele, alt allele")
        #required_ves.add_argument("-m", "--model-h5", type=str, required=True, help="Path model .h5 file")
        #required_ves.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
        #required_ves.add_argument("-op", "--output-prefix", type=str, required=True, help="Output prefix for bigwig files")
   
        #optional_ves.add_argument("-bs", "--batch-size", type=int, default=64, help="batch size to use for prediction")
        #optional_ves.add_argument("-dm","--debug-mode-on", type=int, default=0, help="Use this mode to print the flanks of first five SNP insert locations")
        
        
        # Run TF-Modisco
        
        #motifs_parser._action_groups.pop()
        #required_tfm =  motifs_parser.add_argument_group('required arguments')
        #optional_tfm =  motifs_parser.add_argument_group('optional arguments')
        
        #required_tfm.add_argument("-i", "--h5py", type=str, required=True, help="A legacy h5py file containing the one-hot encoded sequences and shap scores.")
        #required_tfm.add_argument("-n", "--max-seqlets", type=int, required=True, help="The maximum number of seqlets per metacluster.")
        #required_tfm.add_argument("-op", "--output-prefix", type=str, required=True, help="The path to the output file.")

        #optional_tfm.add_argument("-l", "--n-leiden", type=int, default=2, help="The number of Leiden clusterings to perform with different random seeds.")
        #optional_tfm.add_argument("-w", "--window", type=int, default=500, help="The window surrounding the peak center that will be considered for motif discovery.")
        #optional_tfm.add_argument("-v", "--verbose", action="store_true", default=False, help="Controls the amount of output from the code.")

        
        # Pull the arguments
        
        args = parser.parse_args()

        return args

