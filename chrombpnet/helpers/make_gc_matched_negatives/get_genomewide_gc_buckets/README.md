Scripts to find the gc content of genome wide bins.

### Usage

```
chrombpnet_genomewide_gc -g [genome] -o [output_prefix]  -f [inputlen] -s [stride]
```

The above scripts bins the genome into `inputlen` length regions. The bin intervals are spaced at `stride` length.

### Input format:

- genome: Reference geneome fasta.
- output_prefix: Output path prefix to store the output, ".bed" is appended as suffix by the code. If path conatins a directory make sure it already exists. See output format section below.
- inputlen: An integer value for the length of a bin for gc-fraction calculation.
- stride: An integer stride value to space the bins.

### Output format:

The script outputs a tab-seperated file containing chr, start, end and gc_content. These are genome wide regions binned into `inputlen` regions and with a given `stride`.

Genomewide bins for hg38 using source file https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz have been generated
at  http://mitra.stanford.edu/kundaje/anusri/chrombpnet_downloads/genomewide_gc_hg38_stride_1000_inputlen_2114.bed generated with `inputlen` 2114 and `stride` 1000 using this script.

```
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gzip -d GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
chrombpnet_genomewide_gc -g GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -o genomewide_gc_hg38_stride_1000_inputlen_2114 -f 2114 -s 1000
```
