# Scripts to build pwm from bigwig

The following scripts will be used in ChromBPNet to ensure that the bigwigs are shifted correctly.

## Usage

Build a PWM matrix centering at the non-zero entries in a bigwig. 

```
python build_pwm_from_bigwig.py -i [bigwig] -g [genome] -o [output_prefix] -c [chr] -cz [chrom_sizes] -pw [pwm_width]
```
## Example Usage

```
python build_pwm_from_bigwig.py -i unstranded.bw -g GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -o /path/name_of_pwm.png -c chr20 -cz  hg38.chrom.sizes -pw 24
```

## Input

- bigwig: Path to ATAC/DNASE data in bigwig format
- genome: Path to reference genome fasta
- output_prefix: Output prefix path to use for image storage. If prefix includes a directory path make sure it already exists.
- chr: A string value of chromsome to use to build a pwm. This name should be present in both the bigwig file and also should be present in column one of the `chrom_sizes`
- chrom_sizes: Path to a TSV file that has chromosomes in the first column and their sizes in the second column.
- pwm_width: An integer value of PWM width to consider. This defaults to 24.

## Output

Output an PWM image file with the name `output_prefix`.

