# Preprocessing scripts for chrombpnet

The following scripts in this folder are some pre-processing steps to train chrombpnet models

## Requirements

To run these scripts you will need the `samtools` and `bedGraphToBigWig` (from ucsc) tools. These scripts are currently tested **only for bulk ATAC-seq**. If you want to use these for sc-ATAC please make appropriate changes. I will try to add scripts that will work with sc-ATAC soon.

## Bam to Bigwig (+4/-4 shift for ATAC and 0/+1 shift for DNASE)

```
bash bam_to_bigwig.sh [input_bam] [output_prefix] [data_type] [chrom_sizes]
```

The following assumptions are made with this script - make changes accordingly if the assumptions dont hold.

- PE stands for paired-end and SE stands for single-end.
- The DNASE input bams are +1 shifted on negative strand.
- The ATAC input bams are +4 shifted on positive strand and -4 shifted on negative strand.
- When PE is considered we used **filtered bams**. Filtered bams are obtained from the [ENCODE ATAC-seq pipeline][url1]. 
- When SE is considered we start from **unfiltered bams**. We then use samtools flag `780` to do filtering. Refer to the following [link][url2] to understand what this flag means.
- To understand the intuition behind performing this different filtering step on SE and PE please refer to the discussion section [below](#discussion). 
- The script generates an unstranded bigwig - that is forward and reverse strand are not considered seperately but are combined into one.

## Example Usage

```
bash bam_to_bigwig.sh input.bam <output_name> ATAC_PE hg38.chrom.sizes
<example output file> : <output_name>_unstranded.bw
```

## Input Format

- input_bam: must be bam file.
- output_prefix: prefix for output file. The directory in the prefix, if applicable, must already exist.
- data_type: The `data_type` variable takes the following 4 values ATAC_PE, ATAC_SE, DNASE_PE, DNASE_SE.
- chrom_sizes: Tab seperated file with two columns. First column is the chromosome and second column is the chromsome length. (Make sure the bam's use the same reference format as the `chrom_sizes`.)

## Output Format

- Generates a unstranded bigwig file with +4/-4 shift applied on ATAC datasets and 0/+1 shift applied on DNASE datasets.

## Discussion

### Wondering why we are using the  +4/-4 ATAC shift (instead of the popular +4/-5 ATAC shift) and 0/+1 DNASE shift (instead of the popular no shift)? 

Here we show you what a PWM built from unstranded bigwigs look like when considering different shifts. To convert bigwigs to PWM use the scripts provided in the `analysis/` directory. From the below images we see that the signal on forward and revere strand reinforce themselves when we consider +4/-4 shift giving us the Tn5/DNASE-I bias motif PWM that we know of in literature ([Figure 1 in HINT-ATAC][url3] paper).


#### PWM from unstranded ATAC bigwig with no shift

![Image](images/atac_no_shift.png)

#### PWM from  unstranded GM12878 ATAC bigwig with +4/-4 shift on forward/reverse strand respectively 

![Image](images/atac_44_shift.png)

#### PWM from unstranded GM12878 ATAC bigwig with +4/-5 shift on forward/reverse strand respectively

![Image](images/atac_45_shift.png)


#### PWM from unstranded GM12878 DNASE bigwig with no shift 

![Image](images/dnase_no_shift.png)


#### PWM from unstranded GM12878 DNASE bigwig with 0/+1 shift on forward/reverse strand respectively

![Image](images/dnase_01_shift.png)

The images can be reproduced by running the scripts in the `images/directoy`

### What is the intuition behind performing the different filtering step on PE and SE?

When considering paired-end data, duplicates can be marked correctly because we have paired-end reads and hence we can remove them. This is not possible when considering single-end reads and hence we do not remove duplicates (here if we attempt to do this we might be throwing away actual signal information)

[url1]: https://github.com/ENCODE-DCC/atac-seq-pipeline
[url2]: https://broadinstitute.github.io/picard/explain-flags.html
[url3]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2




