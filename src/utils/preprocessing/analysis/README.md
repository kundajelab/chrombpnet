# Scripts to build pwm from bigwig

The following scripts will be used in ChromBPNet to ensure that the bigwigs are shifted correctly.

## Usage

Build a PWM matrix centering at the non-zero entries in a bigwig. 

```
python build_pwm_from_bigwig.py [bigwig] [genome] [output_prefix] [chr] [chr_size] [pwm_width]
```

## Input

- bigwig: Path to ATAC/DNASE data in bigwig format
- genome: Path to reference genome fasta
- output_prefix: Output prefix path to use for image storage. If prefix includes a directory path make sure it already exists.
- chr: A string value of chromsome to use to build a pwm
- chr_size: An integer value of size of the chromsome to be used for building the PWM. The user should make sure that this size wont exceed the chromsome size length.
- pwm_width: An integer value of PWM width to consider. This defaults to 24.

## Output

Output an PWM image file with the name `output_prefix`.