# Bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants

- This repo contains code for the paper [ChromBPNet: Bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants](https://www.biorxiv.org/content/10.1101/2024.12.25.630221v1)
- Wiki with detailed step-by step instructions: https://github.com/kundajelab/chrombpnet/wiki 
- Please contact [Anusri Pampari] (\<first-name\>@stanford.edu) for suggestions and comments. 
- Here is a link to the [slides](https://docs.google.com/presentation/d/1Ow6K8TYN40u7T3ODdo-JRCLuv5fUUacA/edit?usp=sharing&ouid=104820480456877027097&rtpof=true&sd=true), [ISMB talk](https://www.youtube.com/watch?v=3W3JeJvvjLc) and a comprehensive [tutorial](https://github.com/kundajelab/chrombpnet/wiki). Please see the [FAQ](https://github.com/kundajelab/chrombpnet/wiki/FAQ) and file a github [issue](https://github.com/kundajelab/chrombpnet/issues) if you have questions.
- If you are using chrombpnet <= v0.1.3 please refer to the note here - https://github.com/kundajelab/chrombpnet/wiki/Denovo-motif-discovery 
- If you are using chrombpnet repo actively in your project, I strongly recommend adding yourself to the watchers list for updates. Click on the eye symbol (below the star and above the fork symbol to the right). This will keep you informed of all the major updates and bugs posted for this repo.  

Chromatin profiles (DNASE-seq and ATAC-seq) exhibit multi-resolution shapes and spans regulated by co-operative binding of transcription factors (TFs). This complexity is further difficult to mine because of confounding bias from enzymes (DNASE-I/Tn5) used in these assays. Existing methods do not account for this complexity at base-resolution and do not account for enzyme bias correctly, thus missing the high-resolution architecture of these profiles. Here we introduce ChromBPNet to address both these aspects.

ChromBPNet (shown in the image as `Bias-Factorized ChromBPNet`) is a fully convolutional neural network that uses dilated convolutions with residual connections to enable large receptive fields with efficient parameterization. It also performs automatic assay bias correction in two steps, first by learning simple model on chromatin background that captures the enzyme effect (called `Frozen Bias Model` in the image). Then we use this model to regress out the effect of the enzyme from the ATAC-seq/DNASE-seq profiles. This two step process ensures that the sequence component of the ChromBPNet model (called `TF Model`) does not learn enzymatic bias. 

<p align="center">
<img src="images/chrombpnet_arch.png" alt="ChromBPNet" align="center" style="width: 400px;"/>
</p>

> **ChromBPNet 2.x** runs on Keras 3 with the JAX backend instead of TensorFlow. It runs natively
> on CUDA 13 GPUs (H100, B200, RTX PRO 6000 Blackwell), with a single NumPy 2 environment and a native JAX
> DeepSHAP. The commands, inputs and outputs are the same as in 1.x except for the changes listed in the
> [CHANGELOG](CHANGELOG.md), and chrombpnet 1.x models load as they are.

## Table of contents

- [Installation](#installation)
- [GPU notes](#gpu-notes)
- [QuickStart](#quickstart)
- [Optional flags](#optional-flags)
- [Compatibility with chrombpnet 1.x](#compatibility-with-chrombpnet-1x)
- [Validation](#validation)
- [How-to-cite](#how-to-cite)

## Installation

ChromBPNet needs Python >= 3.12 and, for training and interpretation, an NVIDIA GPU. It also runs on the CPU,
including on Apple-silicon Macs, which is fine for testing and small jobs. TensorFlow is not needed. Pick one of
the three ways below.

### 1. pixi (recommended)

[pixi](https://pixi.sh) installs both the Python packages and the command-line tools the pipeline calls
(bedtools, bedGraphToBigWig, samtools, MEME `tomtom`, pango for the PDF reports), all from the committed lock
file:

```
curl -fsSL https://pixi.sh/install.sh | bash
git clone https://github.com/kundajelab/chrombpnet.git
cd chrombpnet
pixi install -e cuda13                    # Linux, NVIDIA driver >= 580
pixi run -e cuda13 gpu-check              # lists the GPU(s) JAX sees
pixi run -e cuda13 chrombpnet pipeline ...
pixi shell -e cuda13                      # or activate the environment and call `chrombpnet` directly
```

| environment | use it for |
|---|---|
| `cuda13` | Linux x86_64 / aarch64 with driver >= 580: H100, B200, RTX PRO 6000 and other Blackwell GPUs |
| `cuda12` | Linux x86_64 with an older driver (>= 525; Blackwell GPUs need >= 570) |
| `default` | CPU only: Linux and macOS arm64 (`pixi install`, `pixi run chrombpnet ...`) |
| `dev`, `cuda13-dev` | the same plus pytest and ruff: `pixi run -e dev test`, `pixi run -e cuda13-dev test-gpu` |

To install a GPU environment on a machine without a GPU (a login node, a CI runner, a container build), prefix
the command with `CONDA_OVERRIDE_CUDA=13.0` (or `12.0`). MEME is only packaged for Linux, so on macOS use
`--tomtom-lite` (see below).

### 2. uv or pip (no conda)

```
git clone https://github.com/kundajelab/chrombpnet.git
cd chrombpnet
uv sync --extra cuda13                    # --extra cuda12 for older drivers, no extra for CPU
uv run chrombpnet --help
# or into an existing virtual environment (Python >= 3.12):
uv pip install -e '.[cuda13]'
```

The chrombpnet releases on PyPI (up to 1.x) are TensorFlow-based. Until 2.x is released there, install it from
source as above. Without conda you have to provide these yourself:

- **bedtools** and **bedGraphToBigWig** (UCSC). `chrombpnet pipeline` and `chrombpnet bias pipeline` need them to
  turn BAM/fragment/tagAlign files into a bigWig, and `chrombpnet prep nonpeaks` needs bedtools.
- **MEME `tomtom`** for motif matching in the reports. Alternatively pass `--tomtom-lite` to use the bundled
  TOMTOM-lite. It is faster and needs no MEME, but reports p-values instead of q-values.
- **pango** for the PDF reports (weasyprint), e.g. `apt install libpango-1.0-0 libpangoft2-1.0-0` or
  `brew install pango`.
- **pyBigWig** has PyPI wheels only for Linux x86_64. Elsewhere it builds from source, which needs a C compiler
  and zlib headers.

### 3. Docker or Apptainer

```
docker run --rm --gpus all ghcr.io/kundajelab/chrombpnet:latest-cuda13 chrombpnet --help
apptainer exec --nv --cleanenv docker://ghcr.io/kundajelab/chrombpnet:latest-cuda13 chrombpnet --help
docker build -t chrombpnet:cuda13 .       # build it yourself; see the Dockerfile header for cuda12 / CPU images
```

The image contains no CUDA toolkit. The host provides only the driver, through the
[NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
(`--gpus`) or `apptainer --nv`. The docker workflow (`.github/workflows/docker.yml`) builds the images for
release tags and pushes them to `ghcr.io/kundajelab/chrombpnet`, with `-cuda13` tags (driver >= 580) and
`-cuda12` tags (older drivers), once the maintainers enable it. Until then, or if there is no image for your
version, build one locally with `docker build`.

Apptainer passes the host environment into the container and binds `$HOME`. The image already ignores the
host's `PYTHONPATH`, `PYTHONHOME` and `~/.local` Python packages. `--cleanenv` also keeps out every other host
variable, including ones you may want. Pass those explicitly, for example
`--env CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"` when a scheduler has picked your GPUs, or
`--env XLA_PYTHON_CLIENT_MEM_FRACTION=0.25`.

## GPU notes

- **Driver.** CUDA 13 needs NVIDIA driver >= 580 (the `nvidia-smi` header shows `CUDA Version: 13.x`). H100
  (sm_90), B200 (sm_100) and RTX PRO 6000 Blackwell (sm_120) are supported natively. With an older driver use the
  `cuda12` environment.
- **No `module load cuda`, no `LD_LIBRARY_PATH`.** JAX brings its own CUDA, cuDNN and NCCL as pip wheels.
  System CUDA libraries found first through `LD_LIBRARY_PATH` shadow them and break start-up.
- **Memory.** Importing chrombpnet sets `XLA_PYTHON_CLIENT_PREALLOCATE=false` when it is unset (the Docker image
  sets it too, and `gpu-check` imports chrombpnet first). JAX then allocates GPU memory as it needs it, instead of
  reserving 75% of the card at start-up. An `export XLA_PYTHON_CLIENT_PREALLOCATE=true` of your own wins.
- **Shared GPUs.** JAX does not return memory it has allocated until the process ends. DeepSHAP picks its batch
  size from the model size and the free GPU memory, so `chrombpnet pipeline` keeps that memory through the
  TF-MoDISco step that follows. On a GPU shared with other jobs, cap the process, e.g.
  `export XLA_PYTHON_CLIENT_MEM_FRACTION=0.25`. The value is a fraction of the card's total memory.
- **Backend.** Importing chrombpnet also sets `KERAS_BACKEND=jax`. If your own code imports `keras` before
  chrombpnet, set `KERAS_BACKEND=jax` in the environment yourself.
- **Clusters.** [`workflows/slurm/`](workflows/slurm/) has a SLURM template for H100/B200 clusters.

## QuickStart

### Bias-factorized ChromBPNet training

The command to train ChromBPNet with pre-trained bias model will look like this:

```
chrombpnet pipeline \
  -ibam /path/to/input.bam \ # only one of ibam, ifrag or itag is accepted
  -ifrag /path/to/input.tsv \ # only one of ibam, ifrag or itag is accepted
  -itag /path/to/input.tagAlign \ # only one of ibam, ifrag or itag is accepted
  -d "ATAC" \
  -g /path/to/hg38.fa \
  -c /path/to/hg38.chrom.sizes \ 
  -p /path/to/peaks.bed \
  -n /path/to/nonpeaks.bed \
  -fl /path/to/fold_0.json \
  -b /path/to/bias.h5 \ 
  -o path/to/output/dir/ \
```

#### Input Format

- `-ibam` or `-ifrag` or `-itag`: input file path with filtered reads in one of bam, fragment or tagalign formats. Example files for supported types - [bam](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_merged.bam), [fragment](https://storage.googleapis.com/chrombpnet_data/input_files/example.fragments.tsv), [tagalign](https://storage.googleapis.com/chrombpnet_data/input_files/example.tagAlign) 
- `-d`: assay type. The following types are supported - "ATAC" or "DNASE"
- `-g`: reference genome fasta file. Example file human reference - [hg38.fa](https://storage.googleapis.com/chrombpnet_data/input_files/hg38.genome.fa)
- `-c`: chromosome and size tab separated file. Example file in human reference - [hg38.chrom.sizes](https://storage.googleapis.com/chrombpnet_data/input_files/hg38.chrom.sizes)
- `-p`: Input peaks in narrowPeak file format, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region 	  is centered at start + summit internally, across all regions. Example file with [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/) dataset - [peaks.bed](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_relaxed_peaks_no_blacklist.bed)
- `-n`: Input nonpeaks (background regions)in narrowPeak file format, and must have 10 columns, with values minimally for chr, start, end and summit 	  	(10th column). Every region is centered at start + summit internally, across all regions. Example file with [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/) dataset - [nonpeaks.bed](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_nonpeaks_no_blacklist.bed). More instructions on how to make your own nonpeak file can be found in the [Preprocessing](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing#generate-non-peaks-background-regions) guide.
- `-fl`: json file showing split of chromosomes for train, test and valid. Example 5 fold jsons for human reference -  [folds](https://zenodo.org/records/7443683/files/folds.zip?download=1) 
- `-b`: Bias model in `.h5` format. Bias models are generally transferable across  assay types following similar protocol. Repository of pre-trained bias models for use [here](https://zenodo.org/records/7443683/files/bias_models.zip?download=1). Instructions to train custom bias model below.
- `-o`: Output directory path

Please find scripts and best practices for preprocssing [here](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing).

#### Output Format

The ouput directory will be populated as follows -

```
models\
	bias_model_scaled.h5
	chrombpnet.h5
	chrombpnet_nobias.h5 (TF-Model i.e model to predict bias corrected accessibility profile) 
logs\
	chrombpnet.log (loss per epoch)
	chrombpnet.log.batch (loss per batch per epoch)
	(..other hyperparameters used in training)
	
auxilary\
	filtered.peaks
	filtered.nonpeaks
	...

evaluation\
	overall_report.pdf
	overall_report.html
	bw_shift_qc.png 
	bias_metrics.json 
	chrombpnet_metrics.json
	chrombpnet_only_peaks.counts_pearsonr.png
	chrombpnet_only_peaks.profile_jsd.png
	chrombpnet_nobias_profile_motifs.pdf
	chrombpnet_nobias_counts_motifs.pdf
	chrombpnet_nobias_max_bias_response.txt
	chrombpnet_nobias.....footprint.png
	...
```
Detailed usage guide with more information on input arguments and the output file formats and how to work with them are provided [here](https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training) and [here](https://github.com/kundajelab/chrombpnet/wiki/Output-format).

For more information, also see:

- [Full documentation list](https://github.com/kundajelab/chrombpnet/wiki)
- [Detailed list of input arguments](https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training)
- [Detailed usage guide with more information on the output file formats and how to work with them](https://github.com/kundajelab/chrombpnet/wiki/Output-format)
- [Best practices for preprocessing](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing)
- [Training tutorial](https://github.com/kundajelab/chrombpnet/wiki/Tutorial)
- [Frequently Asked Questions, FAQ](https://github.com/kundajelab/chrombpnet/wiki/FAQ)
 
## Bias Model training

The command to train a custom bias bias model will look like this:

```
chrombpnet bias pipeline \
  -ibam /path/to/input.bam \ # only one of ibam, ifrag or itag is accepted
  -ifrag /path/to/input.tsv \ # only one of ibam, ifrag or itag is accepted
  -itag /path/to/input.tagAlign \ # only one of ibam, ifrag or itag is accepted
  -d "ATAC" \
  -g /path/to/hg38.fa \
  -c /path/to/hg38.chrom.sizes \ 
  -p /path/to/peaks.bed \
  -n /path/to/nonpeaks.bed \
  -fl /path/to/fold_0.json \
  -b 0.5 \ 
  -o path/to/output/dir/ \
```

#### Input Format

- `-ibam` or `-ifrag` or `-itag`: input file path with filtered reads in one of bam, fragment or tagalign formats. Example files for supported types - [bam](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_merged.bam), [fragment](https://storage.googleapis.com/chrombpnet_data/input_files/example.fragments.tsv), [tagalign](https://storage.googleapis.com/chrombpnet_data/input_files/example.tagAlign) 
- `-d`: assay type.  Following types are supported - "ATAC" or "DNASE"
- `-g`: reference genome fasta file. Example file human reference - [hg38.fa](https://storage.googleapis.com/chrombpnet_data/input_files/hg38.genome.fa)
- `-c`: chromosome and size tab separated file. Example file in human reference - [hg38.chrom.sizes](https://storage.googleapis.com/chrombpnet_data/input_files/hg38.chrom.sizes)
- `-p`: Input peaks in narrowPeak file format, and must have 10 columns, with values minimally for chr, start, end and summit (10th column). Every region 	  is centered at start + summit internally, across all regions. Example file with [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/) dataset - [peaks.bed](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_relaxed_peaks_no_blacklist.bed)
- `-n`: Input nonpeaks (background regions)in narrowPeak file format, and must have 10 columns, with values minimally for chr, start, end and summit 	  	(10th column). Every region is centered at start + summit internally, across all regions. Example file with [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/) dataset - [nonpeaks.bed](https://storage.googleapis.com/chrombpnet_data/input_files/ENCSR868FGK_nonpeaks_no_blacklist.bed)
- `-f`: json file showing split of chromosomes for train, test and valid. Example 5 fold jsons for human reference -  [folds](https://zenodo.org/records/7443683/files/folds.zip?download=1) 
- `-o`: Output directory path

Please find scripts and best practices for preprocessing [here](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing).

#### Output Format

The output directory will be populated as follows -


```
models\
	bias.h5
logs\
	bias.log (loss per epoch)
	bias.log.batch (loss per batch per epoch)
	(..other hyperparameters used in training)
	
intermediates\
	...

evaluation\
        overall_report.html
        overall_report.pdf
	pwm_from_input.png
        k562_epoch_loss.png 
	bias_metrics.json
	bias_only_peaks.counts_pearsonr.png
	bias_only_peaks.profile_jsd.png
	bias_only_nonpeaks.counts_pearsonr.png
	bias_only_nonpeaks.profile_jsd.png
        bias_predictions.h5
	bias_profile.pdf
	bias_counts.pdf
	...
```
Detailed usage guide with more information on the input arguments and output file formats and how to work with them are provided [here](https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training) and [here](https://github.com/kundajelab/chrombpnet/wiki/Output-format).

For more information, also see:

- [Full documentation list](https://github.com/kundajelab/chrombpnet/wiki)
- [Detailed list of input arguments](https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training)
- [Detailed usage guide with more information on the output file formats and how to work with them](https://github.com/kundajelab/chrombpnet/wiki/Output-format)
- [Best practices for preprocessing](https://github.com/kundajelab/chrombpnet/wiki/Preprocessing)
- [Training tutorial](https://github.com/kundajelab/chrombpnet/wiki/Tutorial)
- [Frequently Asked Questions, FAQ](https://github.com/kundajelab/chrombpnet/wiki/FAQ)
 
## Optional flags

These are new in 2.x. The defaults reproduce the 1.x behaviour. `chrombpnet <command> --help` shows which
commands take which flag.

| flag | effect |
|---|---|
| `--optimizer {adam,muon}` | `muon` applies Muon to the dilated convolution kernels and Adam to everything else (experimental) |
| `--muon-lr` | Muon learning rate (default 2e-3) |
| `--ema` | train with an exponential moving average of the weights (momentum 0.999), and evaluate and save the averaged weights |
| `--lr-schedule {constant,cosine}` | learning-rate schedule (default `constant`) |
| `--precision {default,highest,bf16}` | float32 matmul/convolution precision for training: `default` lets the GPU use TF32 as TensorFlow did, `highest` forces full float32, `bf16` trains in mixed bfloat16 |
| `--device {auto,gpu,cpu}` | `gpu` fails immediately if JAX sees no GPU, instead of silently training on the CPU |
| `--interpret-subsample N` | number of peaks used for DeepSHAP and TF-MoDISco in the pipelines (default 30000) |
| `--shap-seed`, `--shap-batch-seqs`, `--shap-precision` | DeepSHAP reference seed (default 1234), sequences per batch (default: automatic), precision (default `auto`: full float32 on CPU, TF32 on GPU as in 1.x) |
| `--modisco-max-seqlets`, `--modisco-window` | TF-MoDISco limits (defaults 50000 and 500) |
| `--tomtom-lite` | match motifs with TOMTOM-lite instead of MEME `tomtom`: much faster, needs no MEME, reports p-values |

## Compatibility with chrombpnet 1.x

- **Old models load directly.** TensorFlow 2.x `.h5` files from chrombpnet 1.x load as they are: `bias.h5`,
  `bias_model_scaled.h5`, `chrombpnet.h5` and `chrombpnet_nobias.h5`, including the
  [pre-trained bias models](https://zenodo.org/records/7443683/files/bias_models.zip?download=1). Use them with
  `-b`, `pred_bw`, `contribs_bw`, `footprints` and so on. The loader replaces the logsumexp `Lambda` layer of
  `chrombpnet.h5` and never runs stored bytecode.
- **Counts-head DeepSHAP needs `chrombpnet_nobias.h5`.** `contribs_bw` refuses the counts head of a full
  `chrombpnet.h5` (1.x or 2.x), because its counts output is a logsumexp of the bias and TF-model heads. Run it
  on `chrombpnet_nobias.h5`, as the pipelines do, or pass `-pc profile`.
- **New models cannot go the other way.** Model files keep their `.h5` names, but Keras 3 writes them.
  TensorFlow 2.x chrombpnet cannot load them, and neither can tools that read the HDF5 weights directly. For
  those, `chrombpnet export --legacy-h5 -m chrombpnet.h5 -o chrombpnet_legacy.h5` writes a copy in the layout of
  the 1.x files, which TF-Keras 2.x, chrombpnet 1.x and bpnet-lite load (see the [CHANGELOG](CHANGELOG.md)).
- **Contribution-score files keep the 1.x layout**: `/raw/seq`, `/shap/seq`, `/projected_shap/seq`, shape
  (N, 4, L), Blosc-compressed. To read them with h5py, `import hdf5plugin` first.
- **Results match 1.x statistically, not bit for bit.** Initial weights and random streams differ from
  TensorFlow for the same seed. DeepSHAP references, the reads used for Tn5/DNase shift estimation and the
  nonpeak subsample behind the outlier thresholds are now seeded. EarlyStopping restores the best epoch. See the
  [CHANGELOG](CHANGELOG.md).

## Validation

[VALIDATION.md](VALIDATION.md) compares 2.x with chrombpnet 1.x on a public ENCODE ATAC-seq dataset (ENCSR763OVZ):
loading ENCODE's released models, predictions and contribution scores against ENCODE's bigWigs, and retraining
with Adam (the default) and Muon. It also describes the test suite and what has not been tested yet.

## How to Cite

If you're using ChromBPNet in your work, please cite as follows:

```
@article {Pampari2024.12.25.630221,
	author = {Pampari, Anusri and Shcherbina, Anna and Kvon, Evgeny and Kosicki, Michael and Nair, Surag and Kundu, Soumya and Kathiria, Arwa S. and Risca, Viviana I. and Kuningas, Kristiina and Alasoo, Kaur and Greenleaf, William James and Pennacchio, Len A. and Kundaje, Anshul},
	title = {ChromBPNet: bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants},
	elocation-id = {2024.12.25.630221},
	year = {2024},
	doi = {10.1101/2024.12.25.630221},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/12/25/2024.12.25.630221},
	eprint = {https://www.biorxiv.org/content/early/2024/12/25/2024.12.25.630221.full.pdf},
	journal = {bioRxiv}
}
```


