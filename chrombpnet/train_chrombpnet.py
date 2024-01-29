import parsers as parsers
import os
from data import DefaultDataFile, get_default_data_path
from data import print_meme_motif_file
import pipelines as pipelines
import copy
import pandas as pd
import logging
logging.getLogger('matplotlib.font_manager').disabled = True

args = parsers.read_parser()

os.makedirs(os.path.join(args.output_dir,"logs"), exist_ok=True)
os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=True)
os.makedirs(os.path.join(args.output_dir,"models"), exist_ok=True)
os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=True)

pipelines.chrombpnet_train_pipeline(args)