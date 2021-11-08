## Preprocessing scripts for chrombpnet

The following scripts in this folder are some pre-processing steps to train chrombpnet models

### Bam to Bigwig (+4/-4 shift for ATAC and 0/+1 shift for DNASE)

```
bash bam_to_bigwig.sh [input_bam] [output_dir] [data_type] [chrom_sizes]
```

The following assumptions are made with this script - make changes accordingly if the assumptions dont hold.

- PE stands for paired-end and SE stands for single-end.
- The DNASE input bams are +1 shifted on negative strand.
- The ATAC input bams are +4 shifted on positive strand and -4 shifted on negative strand.
- When PE is considered we used filtered bams. Filtered bams are obtained from the [ENCODE ATAC-seq pipeline][url1]. 
- When SE is considered we start from unfiltered bams. We then use samtools flag `780` to do filtering. Refer to the following [link][url2] to understand what this flag means.

#### Input Format

- input_bam: must be bam file.
- output_dir: Directory to store the output. The directory must already exist.
- data_type: The `data_type` variable takes the following 4 values ATAC_PE, ATAC_SE, DNASE_PE, DNASE_SE.
- chrom_sizes: Tab seperated file with two columns. First column is the chromosome and second column is the chromsome length. (Make sure the bam's use the same reference format.)

#### Output Format

- Generates a unstranded bigwig file with +4/-4 shift applied on ATAC datasets and 0/+1 shift applied on DNASE datasets.


[url1]: https://github.com/ENCODE-DCC/atac-seq-pipeline
[url2]: https://broadinstitute.github.io/picard/explain-flags.html


### TileDB database generation

This script is needed only if you are using the `tiledb` indexing for training dataset generation. You can skip this step if your using the `batchgen` training dataset generation.


