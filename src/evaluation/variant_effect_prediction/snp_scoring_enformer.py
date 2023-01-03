from snp_generator_seqs import SNPGenerator
from scipy.spatial.distance import jensenshannon
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import pandas as pd
import os
import argparse
import losses
import numpy as np
import pickle as pkl
from tqdm import tqdm
import tensorflow as tf
# Make sure the GPU is enabled 
assert tf.config.list_physical_devices('GPU'), 'Start the colab kernel with GPU: Runtime -> Change runtime type -> GPU'


# Commented out IPython magic to ensure Python compatibility.
import tensorflow_hub as hub
import joblib
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np

# %matplotlib inline
# %config InlineBackend.figure_format = 'retina'


SNP_SCHEMA = ["CHR", "POS0", "REF", "ALT", "META_DATA"]


def fetch_variant_args():
    parser=argparse.ArgumentParser(description="variant effect scoring scripts on SNPS")
    parser.add_argument("-i", "--snp_data", type=str, required=True, help="Path to a tsv output with the following information in columns - chr, position to insert allele (0-based), ref allele, alt allele")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-o","--output_dir", type=str, required=True, help="Path to storing snp effect score predictions from the script, directory should already exist")
    parser.add_argument("-bs","--batch_size", type=int, default=64, help="Batch size to use for model")
    parser.add_argument("-dm","--debug_mode_on", type=int, default=0, help="Use this mode to print the flanks of first five SNP insert locations")
    args = parser.parse_args()
    return args



transform_path = 'gs://dm-enformer/models/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file = '/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa'
clinvar_vcf = 'clinvar.vcf.gz'

# Download targets from Basenji2 dataset 
# Cite: Kelley et al Cross-species regulatory sequence activity prediction. PLoS Comput. Biol. 16, e1008050 (2020).
targets_txt = 'https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt'
df_targets = pd.read_csv(targets_txt, sep='\t')
df_targets.head(3)

"""### Code (double click on the title to show the code)"""

# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)


class EnformerScoreVariantsRaw:

  def __init__(self, tfhub_url, organism='human'):
    self._model = Enformer(tfhub_url)
    self._organism = organism
  
  def predict_on_batch(self, inputs):
    ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
    alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]

    return alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)


class EnformerScoreVariantsNormalized:

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      transform_pipeline = joblib.load(f)
    self._transform = transform_pipeline.steps[0][1]  # StandardScaler.
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)


class EnformerScoreVariantsPCANormalized:

  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      self._transform = joblib.load(f)
    self._num_top_features = num_top_features
    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)[:, :self._num_top_features]


# TODO(avsec): Add feature description: Either PCX, or full names.

# @title `variant_centered_sequences`

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


def variant_generator(vcf_file, gzipped=False):
  """Yields a kipoiseq.dataclasses.Variant for each row in VCF file."""
  def _open(file):
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file)
    
  with _open(vcf_file) as f:
    for line in f:
      if line.startswith('#'):
        continue
      chrom, pos, id, ref, alt_list = line.split('\t')[:5]
      # Split ALT alleles and return individual variants as output.
      for alt in alt_list.split(','):
        yield kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos,
                                           ref=ref, alt=alt, id=id)


def one_hot_encode(sequences):
  batch_data = []
  for sequence in sequences:
      #print(len(sequence))
      batch_data.append(kipoiseq.transforms.functional.one_hot_dna(sequence.upper()).astype(np.float32))
  return np.array(batch_data)


def variant_centered_sequences(vcf_file, sequence_length, gzipped=False,
                               chr_prefix=''):
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))

  for variant in variant_generator(vcf_file, gzipped=gzipped):
    interval = Interval(chr_prefix + variant.chrom,
                        variant.pos, variant.pos)
    interval = interval.resize(sequence_length)
    center = interval.center() - interval.start

    reference = seq_extractor.extract(interval, [], anchor=center)
    alternate = seq_extractor.extract(interval, [variant], anchor=center)

    yield {'inputs': {'ref': one_hot_encode(reference),
                      'alt': one_hot_encode(alternate)},
           'metadata': {'chrom': chr_prefix + variant.chrom,
                        'pos': variant.pos,
                        'id': variant.id,
                        'ref': variant.ref,
                        'alt': variant.alt}}

# @title `plot_tracks`

def plot_tracks(tracks, interval, height=1.5):
  fig, axes = mlp.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
  for ax, (title, y) in zip(axes, tracks.items()):
    ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
    ax.set_title(title)
    #sns.despine(top=True, right=True, bottom=True)
  ax.set_xlabel(str(interval))
  mlp.tight_layout()
  mlp.show()

"""## Make predictions for a genetic sequenece"""


#fasta_extractor = FastaStringExtractor(fasta_file)


# @title Make predictions for an genomic example interval
#target_interval = kipoiseq.Interval('chr11', 35_082_742, 35_083_742)  # @param





def fetch_snp_predictions(snp_regions, inputlen, genome_fasta, batch_size, debug_mode_on=False):
    '''
    Returns model predictions (counts and profile probability predictions) at the given reference and alternate snp alleles.
    Please note that if the SNP location is at the edge - i.e we are unable to form a given inputlen of sequence - we skip predictions at this SNP

    Arguments::
        snp_regions: pandas dataframe with the following columns "CHR", "POS0", "REF", "ALT"
        inputlen: integer representing the input length to use, snp is inserted in the middle
        genome_fasta: path to reference genome
        batch_size: integer value with batch size to use for the model
        debug_mode_on: Takes 0/1 value. Set this  to 1 to print the flanks of first five SNP insert locations. Predictions will be provided only on the these 5 locations.
    
    Returns:
       rsids: Numpy array with (N,) SNP ids. SNP id is a string with the following values "CHR", "POS0", "REF", "ALT" concatenated with delimiter "_". 
            For each of these ids we return the predictions in the lists below. 
       ref_logcount_preds: log count predictions at the reference allele with size (N,)
       alt_logcount_preds: log count predictions at the alternate alele with size (N,)
       ref_prob_preds: profile probability predictions at the reference allele with size (N,outputlen). outputlen depends on the model.
       alt_prob_preds:  profile probability predictions at the alternate allele with size (N,outputlen). outputlen depends on the model.
    '''
    rsids = []
    ref_logcount_preds=[]
    alt_logcount_preds=[]

    # @title Plot tracks
    model = Enformer(model_path)

    # snp sequence generator 
    snp_gen=SNPGenerator(snp_regions=snp_regions,
                        inputlen=inputlen,
                        genome_fasta=genome_fasta,
                        batch_size=batch_size,
                        debug_mode_on=debug_mode_on)

    ref_logcount_preds = np.zeros((len(snp_gen),896))
    alt_logcount_preds = np.zeros((len(snp_gen),896))
    used_ids = []
    locs_list = []

    for i in tqdm(range(len(snp_gen))):
#    for i in range(0,2):

        batch_rsids, variant_locs, ref_seqs, alt_seqs = snp_gen[i]
        if len(ref_seqs) == 0:
            continue
        if len(alt_seqs) == 0:
            continue

        sequence_one_hot = one_hot_encode(ref_seqs)
        predictions = model.predict_on_batch(sequence_one_hot)['human'][0]
        ref_logcount_preds[i,:] = predictions[:, 69]
        del predictions
        del sequence_one_hot

        sequence_one_hot = one_hot_encode(alt_seqs)
        predictions = model.predict_on_batch(sequence_one_hot)['human'][0]       
        alt_logcount_preds[i,:] =  predictions[:, 69]
        del predictions
        del sequence_one_hot

        rsids.extend(batch_rsids)
        used_ids.append(i)
        locs_list.extend(variant_locs)

    return  np.array(used_ids), np.array(rsids), np.array(ref_logcount_preds), np.array(alt_logcount_preds), np.array(locs_list)

def predict_snp_effect_scores(rsids, ref_count_preds, alt_count_preds, ref_prob_preds, alt_prob_preds):
    '''
    Predicts variant effect scores based on model predictions.

    Arguments::
       ref_logcount_preds: log count predictions at the reference allele with size (N,)
       alt_logcount_preds: log count predictions at the alternate alele with size (N,)
       ref_prob_preds: profile probability predictions at the reference allele with size (N,outputlen). outputlen depends on the model.
       alt_prob_preds:  profile probability predictions at the alternate allele with size (N,outputlen). outputlen depends on the model.
    
    Returns:
        log_counts_diff: difference in log count predictions of alternate and reference allele (N,)
        log_probs_diff_abs_sum: Sum of absolute difference in log probability prediction of alternate and reference allele per base. (N,)
        probs_jsd_diff: Jensenshannon distance between probability predictions of alternate and reference allele (N,)
    '''
    log_counts_diff = alt_count_preds - ref_count_preds
    log_probs_diff_abs_sum =  np.sum(np.abs(np.log(alt_prob_preds) -  np.log(ref_prob_preds)),axis=1)*np.sign(log_counts_diff)
    probs_jsd_diff = np.array([jensenshannon(x,y) for x,y in zip(alt_prob_preds, ref_prob_preds)])*np.sign(log_counts_diff)

    return log_counts_diff, log_probs_diff_abs_sum, probs_jsd_diff


if __name__=="__main__":

    args = fetch_variant_args()
    debug_mode_on = args.debug_mode_on


    # load the snp data
    snp_regions=pd.read_csv(args.snp_data,header=None,sep='\t', names=SNP_SCHEMA)
    snp_regions["META_DATA"].fillna('', inplace=True)
    snp_regions['RSID']=snp_regions['CHR'].astype(str)+'_'+snp_regions['POS0'].astype(str)+'_'+snp_regions['REF'].astype(str)+'_'+snp_regions['ALT'].astype('str')+"_"+snp_regions['META_DATA'].astype('str')
    print("printing first 5 rows of the input SNP data provided..")
    print(snp_regions.head(5))

    if debug_mode_on:
        snp_regions = snp_regions.head(5)

    # infer input length
    inputlen=SEQUENCE_LENGTH
    print("input length inferred from the model: ", inputlen)

    # fetch model prediction on snps
    used_ids, rsids, ref_preds, alt_preds, locs_list  = fetch_snp_predictions(snp_regions, inputlen, args.genome, args.batch_size, debug_mode_on)

    # store predictions at snps too - can compute variant effect metrics of your interest - let me know if you find something interesting :)
    data={}
    data["ids"] = used_ids
    data["rsids"] = rsids
    data["locs_list"] = locs_list
    data["ref_logcount_preds"] = ref_preds
    data["alt_logcount_preds"] = alt_preds

    
    pkl.dump(data, open(os.path.join(args.output_dir+"predictions_at_snp.pkl"),'wb'))

    print(ref_preds.shape)
    print(alt_preds.shape)

    ref_prob_preds = ref_preds/ref_preds.sum(axis=1)[:,None]
    alt_prob_preds = alt_preds/alt_preds.sum(axis=1)[:,None]

    print(ref_prob_preds.shape)
    print(alt_prob_preds.shape)

    ref_logcount_preds = np.log(np.squeeze(np.sum(ref_preds, axis = 1))) 
    alt_logcount_preds = np.log(np.squeeze(np.sum(alt_preds, axis = 1)))

    print(ref_logcount_preds.shape)
    print(alt_logcount_preds.shape)


    # find varaint effect scores at snps
    log_counts_diff, log_probs_diff_abs_sum, probs_jsd_diff = predict_snp_effect_scores(rsids, ref_logcount_preds, alt_logcount_preds, ref_prob_preds, alt_prob_preds)

    # unpack rsids to write outputs and write score to output
    snp_effect_scores_pd=pd.DataFrame()
    snp_effect_scores_pd[["CHR", "POS0", "REF", "ALT", "META_DATA"]] = pd.Series(rsids).str.split('_', expand=True)
    snp_effect_scores_pd["log_counts_diff"] = log_counts_diff[used_ids]
    snp_effect_scores_pd["log_probs_diff_abs_sum"] = log_probs_diff_abs_sum[used_ids]
    snp_effect_scores_pd["probs_jsd_diff"] = probs_jsd_diff[used_ids]

    snp_effect_scores_pd.to_csv(os.path.join(args.output_dir, "snp_scores.tsv"), sep="\t", index=False)



