Scripts to find the gc content of genome wide bins.

### Usage

```bash 
python get_genomewide_gc_bins.py -g [genome]  -c [chromosome_sizes_for_reference] -o [output_path_prefix]  -il [inputlen] --s [stride]
```

The above scripts bins the genome into `inputlen` length regions. The bin intervals are spaced at `stride` length.

### Input format:

- genome: Reference geneome fasta.
- chromosome_sizes_for_reference: Tab seperated file with two columns. First column is the chromosome and second column is the chromsome length.
- output_path_prefix: Path and file name to store the output. See output fotmat section below.
- inputlen: An integer value for the length of a bin for gc-fraction calculation.
- stride: An integer stride value to space the bins.

### Output format:

The script outputs a tab-seperated file containing chr, start, end and gc_content. These are genome wide regions binned into `inputlen` regions and with a given `stride`. 
For example look at the following pre-generated file at `/oak/stanford/groups/akundaje/anusri/refs/genomewide_gc_hg38_stride_50_inputlen_1057.bed` generated with `inputlen` 2114 and `stride` 50 using script `/oak/stanford/groups/akundaje/anusri/refs/get_genomewide_gc_bins.sh`.
