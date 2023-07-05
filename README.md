# ChromBPNet-pre-release

This repo is deprecated and can only be used to reproduce results in the linked paper. For the user-facing main page of this repo please refer to - https://github.com/kundajelab/chrombpnet

### Installation

This version works with TensorFlow v2.5.1. Please don't use TF>=v2.6.0 as it has memory leaks that surface late in training! Clone the repo and start by running

`pip install -r requirements.txt` 

in a new conda environment preferably. See [this](https://www.tensorflow.org/install/source#gpu) link to find the appropriate CUDA and cuDNN versions.

### Conventions
- The model that predicts enzyme sequence preference from sequence is called the bias model.
- The model that predicts observed ATAC/DNase-seq signal is called the ChromBPNet model. 
- All bed files are treated as narrowPeak files, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region is centered at start + summit internally, across all files.
- In general there are no constraints or checks on supplied intervals and they are allowed to overlap. When converting output to bigwigs (e.g. model predictions or shap scores), overlapping regions are resolved as follows: the value at a given base is picked from the interval whose summit is closest to the base.
- Effective tn5 shift of +4/-4 is used for ATAC-seq bias models and should be ensured when making bigwigs for training

### Making Base-resolution BigWigs

Base resolution bigwigs for ATAC/DNase-seq represent transposition/digestion events- e.g. given an ATAC-seq fragment from paired end sequencing, there are two transposition (cutting) events, one at each end of the fragment, which would translate to adding +1 to the two coordinates to which the fragment ends map. The bigwigs used in this repo are unstranded, i.e. the `+` and `-` strand ends are combined into one bigwig of "transposition events". 

For ATAC-seq, the effective shift to aligned reads should be +4 on the `+` strand and -4 on the `-` strand. 

### ChromBPNet Model Training

In general, you should download and use a pretrained bias model unless you have a compelling reason to train one yourself. See the `pretrained` directory for bias models.

To train a ChromBPNet model:

```bash
python src/train_chrombpnet.py -g [reference genome fasta] -b [bigwig file] -p [peak bed] -n [non-peak bed] -o [output prefix] -bm [bias model]
```

Peak file should typically contain relaxed peak calls (e.g. overlap peaks from ENCODE pipeline) covering a dynamic range of accessibility. Non-peak file should contain regions that do not overlap with regions in the peak file and typically are GC-matched with peaks (this is not enforced). During training, sequences are jittered by up to 500 bp on either side (controlled by `--max-jitter` parameter).

You should ideally provide at least as many non-peaks as peaks. This is not enforced, but too few non-peaks may break some code. By default, non-peaks are sampled each epoch to contain 1:10 ratio of non-peaks to peaks. This can be modified with the `--negative-sampling-ratio` flag, or if you wish to use all negatives as provided, use the `--no-negative-sampling` flag.

The default model architecture takes in a 2114 bp input, has 512 filters per dilated convolution, 8 dilated convolutions for an output of 1000bp. Due to the model architecture, not all combinations of input size, number of dilated convolutions and output size are acceptable. If you wish to change the input/output sizes, you will have have to calculate a combination that works for your requirement. 

When you run `train_chrombpnet.py`, it also adjusts the bias model’s count outputs to roughly match the sequencing depth of the input sample, and stores the adjusted bias model in the output directory (`output_prefix.adjusted_bias_model.h5`) Use this as the accompanying bias model for downstreams tasks like metrics, prediction, etc.

### Metrics

```bash
python src/metrics.py -g [reference genome fasta] -p [peak bed] -n [non-peak bed] -o [output prefix] -bm [adjusted bias model] -cm [ChromBPNet model] -tc [test chromosome, eg. chr1]
```

Outputs metrics in `output_prefix.metrics.json`.

### Storing Predictions BigWig

```bash
python scripts/predict.py -g [reference genome fasta] -r [regions to predict bed] -o [output prefix] -bm [bias model] -cm [ChromBPNet model] -c [2 column chrom.sizes file] --tqdm 1 [for progress bar] 
```

Runs prediction for regions specified in the regions file and saves 4 outputs: predictions with enzyme bias (which would be most similar to observed track), predictions without bias and one file for each bigwig which contains quantiles of values written (can be used for setting axis limits).

Currently, you need to enter the `scripts` directory to run this.

### Interpretation

```bash
python src/interpret.py -g [reference genome fasta] -r [regions to interpret bed] -m [model (bias or ChromBPNet)] -o [output prefix]
```

This produces the profile and count SHAP scores (both by default, `--profile-or-counts` to choose only one] for each of the specified regions. In addition, it keeps a copy of interpreted regions at `output_prefix.interpreted_regions.bed`. 

You may run this on both the bias and ChromBPNet models.

### Visualizing Interpretation

```bash
python scripts/importance_hdf5_to_bigwig.py -h5 [interpret output h5 (profile or counts)] -r [interpreted regions] -c [2 column chrom.sizes file] -o [bigwig output name] -s [quantiles file output name]
```

This converts the hdf5 output of the interpretation step (counts or profile) into a BigWig, which can be viewed through dynseq tracks on WashU or HiGlass.

Currently, you need to enter the scripts directory to run it. Please submit a PR if you can figure out how to run it without entering the directory...imports are so confusing.

### TF-MoDISco

```bash
python src/run_modisco.py -s [prefix for importance scores path (interpret out prefix)] -p [which set of scores to use (“profile” or “counts”)] -o [output dir]
```

This runs MoDISco on either profile or counts SHAP scores. You can change the maximum number of seqlets used with `--max-seqlets` and the length to which each 
region is center cropped with `--crop`.

A necessary condition for successful ChromBPNet model training is that the MoDISco output on the profile SHAP should not contain any motifs that resemble enzyme bias motifs. 

### Training Bias Models

Don't do it without supervision.

---

To explore more command line options for each command, try `python code_file --help`.

