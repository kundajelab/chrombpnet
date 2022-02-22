Scripts to find the gc content of genome wide bins.

### Usage

```bash 
python get_genomewide_gc_bins.py -g [genome] -o [output_path]  -f [inputlen] -s [stride]
```

The above scripts bins the genome into `inputlen` length regions. The bin intervals are spaced at `stride` length.

### Input format:

- genome: Reference geneome fasta.
- output_path: Path and file name to store the output. See output format section below.
- inputlen: An integer value for the length of a bin for gc-fraction calculation.
- stride: An integer stride value to space the bins.

### Output format:

The script outputs a tab-seperated file containing chr, start, end and gc_content. These are genome wide regions binned into `inputlen` regions and with a given `stride`. 
For example look at the following pre-generated file at http://mitra.stanford.edu/kundaje/anusri/chrombpnet_downloads/genomewide_gc_hg38_stride_1000_inputlen_2114.bed generated with `inputlen` 2114 and `stride` 1000 using this script.