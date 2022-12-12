# Bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants

- This repo contains code for the `Bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants` paper Anusri Pampari*, Anna Shcherbina*, Anshul Kundaje. (*authors contributed equally to the work)  
- Please contact [Anusri Pampari] (\<first-name\>@stanford.edu) for suggestions and comments. 

Chromatin profiles (DNASE-seq and ATAC-seq) exhibit multi-resolution shapes and spans regulated by co-operative binding of transcription factors (TFs). This complexity is further difficult to mine because of confounding bias from enzymes (DNASE-I/Tn5) used in these assays. Existing methods do not account for this complexity at base-resolution and do not account for enzyme bias correctly, thus missing the high-resolution architecture of these profile. Here we introduce ChromBPNet to address both these aspects.

ChromBPNet (shown in the image as `chrombpnet model`) is a fully convolutional neural network that uses dilated convolutions with residual connections to enable large receptive fields with efficient parameterization. It also performs automatic assay bias correction in two steps, first by learning simple model on chromatin background that captures the enzyme effect (called `bias model` in the image). Then we use this model to regress out the effect of the enzyme from the ATAC-seq/DNASE-seq profiles. This two step process ensures that the sequence component of the ChromBPNet model (called `sequence model`) does not learn enzymatic bias. 

![Image](images/chrombpnet_arch.pdf)

## Table of contents

- [Installation](#installation)
- [QuickStart](#quickstart)
- [How-to-cite](#how-to-cite)

## Installation

This section will discuss the packages needed to train a ChromBPNet model. Firstly, it is recommended that you use a GPU for model training and have the necessary NVIDIA drivers and CUDA already installed. You can verify that your machine is set up to use GPU's properly by executing the `nvidia-smi` command and ensuring that the command returns information about your system GPU(s) (rather than an error). Secondly there are two ways to ensure you have the necessary packages to train ChromBPNet  models which we detail below,

### 1. Installation setup through Docker

Download and install the latest version of Docker for your platform. Here is the link for the installers -<a href="https://docs.docker.com/get-docker/">Docker Installers</a>.  Run the docker run command below to open a environment with all the packages installed and do `cd chrombpnet` to start running the tutorial.

> **Note:**
> To access your system GPU's from within the docker container, you must have [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) installed on your host machine.

```
docker run -it --rm --memory=100g --gpus device=0  kundajelab/chrombpnet:dev
```

### 2. Installation setup throuh github within a conda environment 

Create a clean conda environment with python >=3.8 
```
conda create -n chrombpnet python=3.8
conda activate chrombpnet
```

Install non-Python  requirements via conda
```
conda install -y -c conda-forge -c bioconda samtools bedtools ucsc-bedgraphtobigwig pybigwig meme
```

Git clone the staging branch of chrombpnet and install via pip

```
git clone https://github.com/kundajelab/chrombpnet.git
pip install -e chrombpnet
```
	
## QuickStart


The command to train ChromBPNet will look like this:

```
train_chrombpnet_model.sh \
  -i ../tests/data/merged_unsorted.bam \
  -d ATAC \
  -g ../tests/data/hg38.fa \
  -p ../tests/data/overlap.bed.gz \
  -c ../tests/data/hg38.chrom.sizes \ 
  -n ../tests/outputs/data/negatives_with_summit.bed \
  -f ../tests/outputs/data/splits/fold_0.json \
  -b ../tests/outputs/models/bias_model/bias.h5 \ 
  -o testing/ \
```

### Input Format


### Output Format


For a optimal bias correction, please refer to our detailed tutorial and FAQs at our [wiki page](https://github.com/kundajelab/chrombpnet/wiki). 
For other downstream tools refer to Wiki.

## How to Cite
