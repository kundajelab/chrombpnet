
# Scripts to do MODSICO on deepshap output of ChromBPNet

The scripts in this folder do the following three steps (1)  Do de-novo motif discovery on the deepshap output of chrombpnet using MODISCO (2) Annotate the motifs using TOMTOM and (3) Summarize the output to a html format.

## Usage

```
bash run.sh [scores_prefix] [output_dir] [score_type] [seqlets] [crop] [meme_db] [meme_logos] [vier_logos] [vier_html] [html_link]
```

The following assumptions are made with this script - make changes accordingly if the assumptions dont hold.

- The following scripts are used on the output of `/src/evaluation/interpret/interpret.py`. 

## Example Usage

```
bash run.sh /path/to/deepshap_scores/ /path/to/store/output/ counts/profiles 200000 1000 /oak/stanford/groups/akundaje/soumyak/motifs/motifs.meme.txt /oak/stanford/groups/akundaje/soumyak/motifs/pfms/ /oak/stanford/groups/akundaje/projects/chromatin-atlas/vierstra_logos/ http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/vierstra_logos/  http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper/testing/
```

## Input Format

- scores_prefix: This is the `output_prefix` used with `/src/evaluation/interpret/interpret.py`. 
- score_type: This is either set to `counts` or `profile`.
- output_dir: Path to a directory to store the output files. The script assumes that the directory already exists. Look at the output format section below to understand the files generated.
- seqlets: Number of seqlets to use for modisco run. If using the most recent dev version of MODISCO this can be set to 200K. If using older version set to 50K. You can test the working of the script with a much smaller value - as this decides the runtime of the script.
- crop: An integer value representing the crop length to use on the chrombpnet input. In chrombpnet we get contribution scores for 2114 length input but we will run modisco on only 1000 length input. So the default value for this parameter is set to 1000. We do this to avoid catching "AT" rich nucleosome motifs that occur more frequently on the flanks of the 2114 length input.
- meme_db: A txt file containing the meme motifs letter probability matrix. Example file - `/oak/stanford/groups/akundaje/soumyak/motifs/motifs.meme.txt`
- meme_logos: A directory containing meme pfms - Example file - `/oak/stanford/groups/akundaje/soumyak/motifs/pfms`
- vier_logos: A directory with images for the pfms provided in `meme_logos`. If the image is not already present the scripts creates that image.
- vier_html: An html link to a folder to host the images in `vier_logos`. These links are used in `score_type.motifs.html` page created below in output format section.
- html_link: An html prefix link that is hostable. The folder `untrimmed_logos` should be hostable with this html as prefix.

## Output Format

- modisco_results_allChroms_`score_type`.hdf5: Modisco output `hdf5` file.
- seqlets_`score_type`.txt: 
- `score_type`.tomtom.tsv: Tomom output.
- `score_type`.motifs.html: An html link where the modisco motifs are ranked based on their frequency. This html link will host html images provided in `html_link` and `vier_html`. 
- `untrimmed_logos`: This directory will be created in `output_dir` path provided. This directory will store the motif images that are not trimmed.
- `trimmed_logos`: This directory will be created in `output_dir` path provided. This directory wil store the motif images that are trimmed based on the default threshold of 0.3 used in the scripts. 
