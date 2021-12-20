# ChromBPNet: Deep learning models of base-resolution chromatin profiles
                This repo is under construction. Please check back in a week

- This repo contains code for the chrombpnet paper Anusri Pampari*, Anna Shcherbina*, Anshul Kundaje. (*authors contributed equally to the work)  
- General queries/thoughts have been addressed in the discussion section below.
- Please contact [Anusri Pampari] (\<first-name\>@stanford.edu)  for suggestions and comments. More instructions about reporting bugs detailed below.
- Authors would like to thank Avanti Shrikumar and Surag Nair for their help with the project.

## Quick Links

- [About](#chrombpnet)
- [Installation](#installation)
- [Tutorial](#tutorial-on-how-to-train-chrombpnet-models)
    - [Preprocessing](#preprocessing)
    - [Train and Evaluate Bias Model](#train-and-evaluate-bias-model)
    - [Train and Evaluate ChromBPNet Model](#train-and-evaluate-chrombpnet-model)


##  ChromBPNet

Chromatin profiles (DNASE-seq and ATAC-seq) exhibit multi-resolution shapes and spans regulated by co-operative binding of transcription factors (TFs). This complexity is further difficult to mine because of confounding bias from enzymes (DNASE-I/Tn5) used in these assays. Existing methods do not account for this complexity at base-resolution and do not account for enzyme bias correctly, thus missing the high-resolution architecture of these profile. Here we introduce ChromBPNet to address both these aspects.

ChromBPNet (shown in the image as `chrombpnet model`) is a fully convolutional neural network that uses dilated convolutions with residual connections to enable large receptive fields with efficient parameterization. It also performs automatic assay bias correction in two steps, first by learning simple model on chromatin background that captures the enzyme effect (called `bias model` in the image). Then we use this model to regress out the effect of the enzyme from the ATAC-seq/Dnase-seq profiles. This two step process ensures that the sequence component of the ChromBPNet model (called `sequence model`) does not learn enzymatic bias. 

![Image](images/chrombpnet_arch.PNG)

If you are interested in learning more about the detailed architectures used, please refer to the following architecture files - 

- bias model: https://github.com/kundajelab/chrombpnet/blob/master/src/training/models/bpnet_model.py
- chrombpnet model: https://github.com/kundajelab/chrombpnet/blob/master/src/training/models/chrombpnetwith_bias_model.py.

## Installation

This section will discuss the packages needed to train a chrombpnet model. Firstly, it is recommended that you use a GPU and have the necessary NVIDIA drivers already installed and setup as chrombpnet model training is faster on a GPU. Secondly there are two ways to ensure you have the necessary packages to run train chrombpnet models which we detail below,

### 1. Setup through docker


### 2. Setup through Conda

```
pip install -r requirements.txt
```
in a new conda environment preferably. See this link to find the appropriate CUDA and cuDNN versions

TODO -
setup conda environment
setup docker
test setup in a new environment
	
##  Tutorial on how to train chrombpnet models

Here we provide a step-by-step guide to training and evaluating chrombpnet models using the GM12878 ATAC-seq data (ENCSR095QNB) [here][url1] whic is a is bulk ATAC-seq data.

###  Preprocessing

#### Step 1: Download experimental data

We will first start by creating a directory named `data/downloads` and downloading the corresponding files (bams and peak files) for ENCSR095QNB ENCODE dataset using the commands in the bash script below. 

```
mkdir data
mkdir data/downloads
bash step1_download_bams_and_peaks.sh data/downloads
```

Following are some things to keep in mind when using custom datasets/downloads -
- For bulk ATAC-seq/DNASE-seq dataset we use the latest ENCODE ATAC-seq protocol https://github.com/ENCODE-DCC/atac-seq-pipeline. The pipeline outputs both stringent (IDR peaks) and relaxed (Overlap peak) thresholding of peaks across replicates and we use the relaxed thresholding of peaks.
- If you are downloading the data from the ENCODE portal you can download the peaks flagged default for ATAC-seq datasets. For DNASE-seq datasets you might have to use the MACS2 protocol to call peaks on the filtered bams (TODO - provide scripts)
- For paired end data we download the filtered bams output from the pipeline and for single-end data we download the unfiltered bams from the pipeline. Please refer to the documentation below to understand the reason for this difference `src/helpers/preprocessing/`
- TODO - add notes on how will this be different for scATAC


#### Step 2: Make Bigwigs from bam files (IMPORTANT STEP! PLEASE READ CAREFULLY)

We will now create unstranded bigwigs (i.e. the + and - strand ends are combined into one bigwig) from the downloaded bam files (from step1) using the command below. This script uses the following two commands (1) `src/helpers/preprocessing/bam_to_bigwig.sh`: Considers that the given bam files are *unshifted* and does a shift of +4/-4 and (2) `src/helpers/preprocessing/analysis/build_pwm_from_bigwig.py` generates an image of the bias motif and does a sanity check if the shift is correct/incorrect. 

```
bash step2_make_bigwigs_from_bams.sh data/downloads/merged.bam data/downloads/ ATAC_PE data/downloads/hg38.fa data/downloads/hg38.chrom.sizes
```

After running this command open the `bias_pwm.png` image generated by the script in `data/downloads` folder. You will see the following tn5 motif PWM for this dataset.

![Image](images/bias_pwm.png)

Following are some things to keep in mind when using custom datasets/downloads -

- **IMPORTANT NOTE 1:** If you are running these commands on custom experimental bam - *read the documentation* in the directory `$src/helpers/preprocessing/` to make sure you are using the script correctly. Next use this script `$src/helpers/preprocessing/analysis/` and make sure that the script throws no warnings and that you see Tn5 or DNase-I bias pwm in `bias_pwm.png`. If you do not see this it is likely that the bam's that you provided have a shift of some kind. Please provide only unshifted bams to the script. **Do not proceed further if you do not see a Tn5 or DNase-I motif after this step.** 

- **IMPORTANT NOTE 2:** If you are running the pipeline on custom generated bigwigs (without using `$src/helpers/preprocessing/`) make sure the bigwigs are unstranded and make sure the shifts are done correctly. To check this `$src/helpers/preprocessing/analysis/` run the scripts in this directory and make sure you see the Tn5 or DNase-I bias pwm in `bias_pwm.png`. **Do not proceed further if you do not see a Tn5 or DNase-I motif after this step.** 

#### Step 3: Generate background regions gc-matched with the peaks

Here we will generate non-peak background regions that GC-match with the peak regions. We will use the regions to train and evaluate a bias model. We will also use these regions in model training and as background regions to get marginal footprints.

First we will start by dividing the entire genome into overlapping bins of `inputlen` regions. ChromBpnet models are trained on `inputlen` of 2114. This step takes several hours (Sorry never got to speeding this step since its a one-time step, if you have ways to speed this step please contribute.) For conveinve the human genome (hg38) data we use is provided here -  /oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_2114.bed. But if you want to make this file from scratch follow these steps

```
python src/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets/get_genomewide_gc_bins.py -g data/hg38.fa -c data/hg38.chrom.sizes -o data/genomewide_gc_hg38_stride_50_inputlen_2114.bed -il 2114 --s 50
```

Now we will find the regions from this genome-wide bucket that do not fall in overlap-peaks or blacklist bed but have same gc-distribution as the overlap-peak set

```
bash step3_get_background_regions.sh data/hg38.fa data/hg38.chrom.sizes data/blacklist.bed.gz data/overlap.bed.gz 2114 data/ /oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_2114.bed
```

TODO - provide link for all the documentation

###  Train and Evaluate Bias Model

We are now ready to train a bias model! We will first define how we want to split the dataset by menitoning out splits in `splits.py` and generating the corresponding json. If you want to mention custom splits please edit the `splits.py` file directly. This will generate five json files each with a different fold information. You can choose any fold for model training further. For the purpose of this tutorial we will use `fold_0.json`

```
mkdir data/splits
python src/helpers/make_chr_splits/splits.py -o data/splits
```

#### Step 4: Train Bias Model

This creates a directory `output/bias_model` and stores output here
This step consists of three steps - (1) Generate hyperparmeters file for bias model (2) Train bias model and (3) Predict and evaluate bias model. It is very improtant to read this section about the bias hyper-paramters - explain the threshold factor

```
bash step4_train_bias_model.sh data/shifted.sorted.bam.chrombpnet.unstranded.bw data/overlap.bed.gz data/negatives_with_summit.bed data/hg38.fa 0.9 2114 1000 data/splits/fold_0.json

```

How to interpret the output information. Important information about how bias model can catch motifs. So we will do Step 5 - If bias model catches motifs come back here and reduce the threshold and retrain. Careful not to use too few peaks. What worked for me best? 

#### Step 5: Interpret bias model

We can do this on the entire peak set - but we will do this on some subsampled data for a validation step.
subsample some peaks for interpretation and modisco - this is mostly for time purposes - we will subsample to 30K peaks. We will make sure that the peaks we choose dont overlap with the blacklist regions. 

```
inputlen=2114
blacklist_region=data/blacklist.bed.gz
chrom_sizes=data/hg38.chrom.sizes
overlap_peak=data/overlap.bed.gz

mkdir data/subsample_overlap

flank_size=$(( inputlen/2 ))
bedtools slop -i $blacklist_region -g $chrom_sizes -b $flank_size > data/subsample_overlap/temp.txt
bedtools intersect -v -a $overlap_peak -b data/subsample_overlap/temp.txt | shuf  > data/subsample_overlap/temp_n.txt
shuf -n 30000 data/subsample_overlap/temp_n.txt > data/subsample_overlap/30K.subsample.overlap.bed
rm  data/subsample_overlap/temp.txt
rm data/subsample_overlap/temp_n.txt
```

Interpretation 

```
bash step5_interpret_bias_model.sh  data/hg38.fa data/subsample_overlap/30K.subsample.overlap.bed  output/bias_model/model.0.h5
```

This will generate images - make sure you see only tn5 motif and no bias model.

#### Step 6: Scale bias model (IMPORTANT IF THE BIAS MODEL IS BEING TRANSFERED FROM A DIFFERENT EXPERIMENT)

This step is needed if we are transferring a bias model - this will normalize for the difference in read depth If the bias model is trained on an experiment on one read depth versus being used on another read depth. We are testing how well bias models generally transfer - one quick check is build the bias PWM on both the orginal experiment and transfer experiments and make sure the tn5 motifs looks similar.


similar to bias odel training above additionally we will priovide the bias model that we want to train as an argument.
) Generate hyperparmeters file for bias model (2) Train bias model and

For the purpose of this tutorial we will performing scaling of our model in step 6 to show a demo. In reality this si not needed unless you are transferring a bias model.

```
bash step4_train_bias_model.sh data/shifted.sorted.bam.chrombpnet.unstranded.bw data/overlap.bed.gz data/negatives_with_summit.bed data/hg38.fa 0.9 2114 1000 data/splits/fold_0.json output/bias_model/model.0.h5

```

TODO - 
make output dir and provide as input for all the above
Test this step add scale model in models

###  Train and Evaluate ChromBPNet Model

Now that we have a bias model we will use it to regress out the bias from the ATAC-seq and DNAE-seq signal so that the sequence component of the model can have bias free signal.



#### Step 7: Train ChromBPNet Model (This step will also generate the sequence model)


```
bash step7_train_chrombpnet_model.sh data/shifted.sorted.bam.chrombpnet.unstranded.bw data/overlap.bed.gz data/negatives_with_summit.bed data/hg38.fa 0.9 2114 1000 data/splits/fold_0.json output/bias_model/model.0.h5

```

This will output 2 models one is chrombpnet model (this is the model with bias) and the other is the sequence model (this is the moel without the bias)

#### Step 8: Interpret chrombpnet model and sequence model

We will now interpret on the same peak regions that we used for interpreting the bias model.

```
bash step5_interpret_bias_model.sh  data/hg38.fa data/subsample_overlap/30K.subsample.overlap.bed  output/chrombpnet_model/model_wo_bias.0.h5
```

In the images produced here you should not see the bias motif and should only see the cell-type specfic motifs in the profile modisco.

#### Step 9: Marginal footprinting using bias, chrombpnet and sequence model

This is another quick test to make sure the bias correction is working - do footprinting on the following motifs

Run this on bias model, chrombpnet model,
```

```

