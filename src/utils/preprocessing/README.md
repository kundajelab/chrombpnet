# Preprocessing scripts for chrombpnet

The following scripts in this folder are some pre-processing steps to train chrombpnet models

## Bam to Bigwig (+4/-4 shift for ATAC and 0/+1 shift for DNASE)

```
bash bam_to_bigwig.sh [input_bam] [output_dir] [data_type] [chrom_sizes]
```

The following assumptions are made with this script - make changes accordingly if the assumptions dont hold.

- PE stands for paired-end and SE stands for single-end.
- The DNASE input bams are +1 shifted on negative strand.
- The ATAC input bams are +4 shifted on positive strand and -4 shifted on negative strand.
- When PE is considered we used filtered bams. Filtered bams are obtained from the [ENCODE ATAC-seq pipeline][url1]. 
- When SE is considered we start from unfiltered bams. We then use samtools flag `780` to do filtering. Refer to the following [link][url2] to understand what this flag means.
- Scripts are assumed to be run from the main `chrombpnet_paper` repo. You can run these scripts from within this directory too - just update the python paths in the script.

## Input Format

- input_bam: must be bam file.
- output_dir: Directory to store the output. The directory must already exist.
- data_type: The `data_type` variable takes the following 4 values ATAC_PE, ATAC_SE, DNASE_PE, DNASE_SE.
- chrom_sizes: Tab seperated file with two columns. First column is the chromosome and second column is the chromsome length. (Make sure the bam's use the same reference format.)

## Output Format

- Generates a unstranded bigwig file with +4/-4 shift applied on ATAC datasets and 0/+1 shift applied on DNASE datasets.

## Discussion

### Wondering why we are using the  +4/-4 ATAC shift (instead of the popular +4/-5 ATAC shift) and 0/+1 DNASE shift (instead of the popular no shift)? 

Here we show you what a PWM built from unstranded bigwigs look like when considering different shifts. To convert bigwigs to PWM use the scripts provided in the `analysis/` directory. From the below images we see that the signal on forward and revere strand reinforce themselves when we consider +4/-4 shift giving us the Tn5/DNASE-I shift that we know of in literature (Figure 1 in HINT-ATAC [url3] papar)

##### PWM from unstranded ATAC bigwig with no shift

![Image](images/atac_no_shift.png)

##### PWM from  unstranded GM12878 ATAC bigwig with +4/-4 shift on forward/reverse strand respectively 

![Image](images/atac_44_shift.png)

##### PWM from unstranded GM12878 ATAC bigwig with +4/-5 shift on forward/reverse strand respectively

![Image](images/atac_45_shift.png)


##### PWM from unstranded GM12878 DNASE bigwig with no shift 

![Image](images/dnase_no_shift.png)


##### PWM from unstranded GM12878 DNASE bigwig with 0/+1 shift on forward/reverse strand respectively

![Image](images/dnase_01_shift.png)

[url1]: https://github.com/ENCODE-DCC/atac-seq-pipeline
[url2]: https://broadinstitute.github.io/picard/explain-flags.html
[url3]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2


### TileDB database generation

This script is needed only if you are using the `tiledb` indexing for training dataset generation. You can skip this step if your using the `batchgen` training dataset generation.


