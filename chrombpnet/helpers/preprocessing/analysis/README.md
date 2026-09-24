# Scripts to build pwm from bigwig

The following scripts will be used in ChromBPNet to ensure that the bigwigs are shifted correctly (`chrombpnet pipeline` and `chrombpnet bias pipeline` write the plot to `evaluation/bw_shift_qc.png`).

## Usage

Build a PWM matrix centering at the non-zero entries in a bigwig. 

```
python -m chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig -i [bigwig] -g [genome] -op [output_prefix] -cr [chr] -c [chrom_sizes] -pw [pwm_width]
```
## Example Usage

```
python -m chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig -i unstranded.bw -g GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -op /path/name_of_pwm -cr chr20 -c hg38.chrom.sizes -pw 24
```

## Input

- bigwig: Path to ATAC/DNASE data in bigwig format
- genome: Path to reference genome fasta
- output_prefix: Output prefix of file name to use for image storage. If prefix includes a directory path make sure it already exists. Code will append .png suffix to provided argument.
- chr: A string value of chromsome to use to build a pwm. This name should be present in both the bigwig file and also should be present in column one of the `chrom_sizes`
- chrom_sizes: Path to a TSV file that has chromosomes in the first column and their sizes in the second column.
- pwm_width: An integer value of PWM width to consider. This defaults to 24.

## Output

Output an PWM image file with the name `output_prefix.png`. The logo is the information-content scaled PWM (`chrombpnet.utils.viz_sequence`, the TF-MoDISco 0.5 logo code vendored into chrombpnet).
