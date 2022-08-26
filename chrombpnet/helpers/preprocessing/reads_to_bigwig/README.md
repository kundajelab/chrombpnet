# Preprocessing scripts for chrombpnet

The scripts in this folder are pre-processing steps to convert input reads to bigwig format for trainining chrombpnet models.

## Requirements

To run these scripts you will need the `samtools` and `bedGraphToBigWig` (from ucsc) tools. These scripts are tested for bulk and single-cell ATAC-seq, and for bulk DNase-seq.

## BAM/fragment file/tagAlign file to Bigwig

We convert input BAMS to appropriately shifted (+4/-4 shift for ATAC and 0/+1 shift for DNASE) Bigwigs consistent with our training pipeline.

```
usage: reads_to_bigwig.py [-h] -g GENOME (-ibam INPUT_BAM_FILE | -ifrag INPUT_FRAGMENT_FILE | -itag INPUT_TAGALIGN_FILE) -c CHROM_SIZES -o OUTPUT_PREFIX -d {ATAC,DNASE}
                          [-p PLUS_SHIFT] [-m MINUS_SHIFT] [--ATAC-ref-path ATAC_REF_PATH] [--DNASE-ref-path DNASE_REF_PATH] [--num-samples NUM_SAMPLES]

Automatically detect enzyme shift of input BAM/fragment/tagAlign File

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        reference genome file
  -ibam INPUT_BAM_FILE, --input-bam-file INPUT_BAM_FILE
                        Input BAM file
  -ifrag INPUT_FRAGMENT_FILE, --input-fragment-file INPUT_FRAGMENT_FILE
                        Input fragment file
  -itag INPUT_TAGALIGN_FILE, --input-tagalign-file INPUT_TAGALIGN_FILE
                        Input tagAlign file
  -c CHROM_SIZES, --chrom-sizes CHROM_SIZES
                        Chrom sizes file
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix (path/to/prefix)
  -d {ATAC,DNASE}, --data-type {ATAC,DNASE}
                        assay type
  -p PLUS_SHIFT, --plus-shift PLUS_SHIFT
                        Plus strand shift applied to reads. Estimated if not specified
  -m MINUS_SHIFT, --minus-shift MINUS_SHIFT
                        Minus strand shift applied to reads. Estimated if not specified
```

Please supply one of BAM(`-ibam`)/fragment file(`-ifrag`)/tagAlign file(`-itag`) as input. The script generates an unstranded bigwig- forward and reverse strands are combined with appropriate shifting (+4/-4 for ATAC and 0/+1 for DNase). Output is stored at `{OUTPUT_PREFIX}_unstranded.bw`. The directory in the prefix, if applicable, must already exist.

If supplying a fragment file, it should minimally have 3 columns for chr, start and end. Each line must represent a fragment with Tn5 transposition events at both ends.

If supplying a fragment file, it should minimally contain 3 columns for chr, start and end, and a 6th column with the strand.

The CHROM_SIZES file should be a tab-separated file with two columns. First column is the chromosome and second column is the chromsome length. Make sure the bam's use the same reference format as the `CHROM_SIZES`.

### Automatic shift detection

Most ATAC-seq (single-cell and bulk) pipelines shift Tn5 reads by +4/-5 by default. However when combining analyses with other tools, the effective shift can be different than +4/-5 due to off-by-one errors. Our script handles such cases by default by automatically detecting the enzyme shift (for ATAC and DNase) and correcting it appropriately (+4/-4 for ATAC and 0/+1 for DNase) for consistency with the training pipeline.

In rare cases, you may see an error such as "Input file shifts inconsistent" or "Input shift is non-standard". In such cases, if you know the actual shift for your input file (typically +0/+0 for BAMs, and +4/-5 for ATAC fragment/tagAligns) you can supply them using the `--plus-shift` and `--minus-shift` flags. However, if you are uncertain, please reach out to us by submitting an [Issue](https://github.com/kundajelab/chrombpnet/issues).

## Discussion

### Wondering why we are using the  +4/-4 ATAC shift (instead of the popular +4/-5 ATAC shift) and 0/+1 DNASE shift (instead of the popular no shift)? 

**TODO**: Anusri now that we are not using PWMs, maybe provide a different argument?
Here we show you what a PWM built from unstranded bigwigs look like when considering different shifts. To convert bigwigs to PWM use the scripts provided in the `analysis/` directory. From the below images we see that the signal on forward and revere strand reinforce themselves when we consider +4/-4 shift giving us the Tn5/DNASE-I bias motif PWM that we know of in literature ([Figure 1 in HINT-ATAC][url1] paper).

### Suggestions for preparing input data

Please ensure your input files are filtered for quality metrics and duplicates are removed appropriately. For single-end data, we recommend not removing duplicates. Fragment files obtained from Chromap, cellranger and similar tools can be supplied as is.

[url1]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2




