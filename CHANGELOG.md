# Changelog

##  Version - 1.5
- Fixed issue #150, regions_used not found while generating bigwigs from impotance h5s

## Version - 1.4
- (MAJOR) Bug in chrombpnet modisco_motifs command. seqlets was limited to 50000. If users wanted to change it to 1 million this did not happen.
- Filter peaks at edges for pred_bw command and bias pipleline. So bias evaluation now done on these filtered peaks.
- Preprocessing deafulted to use unix sort. Provided option to switch to bedtools sort.
- Provided option to use filter chromosomes option in preprocessing.

## Version - 1.3 - Inworks - 2022-12-11
- Changed pipelines to use modisco-lite, old modisco will soon be removed
- Added automatic shift scripts to repo and integrated with the pipeline
- modisco now outputs both html and pdf. pdf can be shared with anyone.
- Use ATAC, DNASE not ATAC_PE, DNASE_SE anymore.
- Simplyfing workflows to include only two main workflows, chrombpnet_train_tf_model, chrombpnet_train_bias_model
- Restructuring README, moving tutorial to additional documentation and introducing FAQ, and only two pipelines (chrombpnet_train_tf_model, chrombpnet_train_bias_model)

## [Unreleased] - 2022-02-28
- Typo fix - (PR#31-36)
- We can now specify ylimit in marginal footprinting plots (PR#27)
- PR#30 merged  to do fast genome wide gc binning and bug fix to ensure case-sensitve GC calculation. Adds a unit test to make checks.

## [Unreleased] - 2022-02-11
- Updated modisco version - (PR#25)

## [Unreleased] - 2022-02-07
- Typo fix - (PR#23-24)

## [Unreleased] - 2022-01-31
- Marginal footprinting no longer hardcoded to take only chr1, it now inputs the fold json. (PR#20)
- Invivo footprinting is made compatible with the new input repo changes. (PR#22)
- Pandas read bugs fixed for gc-matching scripts. (PR#15)
- Tutorial now exits when one part of the pipeline breaks instead of proceeding further. (PR#18)
- Added the fast binning scripts at src/helpers/make_gc_matched_negatives/get_genomewide_gc_buckets (PR#11)
- Added pseudocount in metrics.py (PR#13) - needs further review.

## [Unreleased] - 2022-01-24
- Seed setting fixed for training. Seed was set improperly - preventing shuffling for every epoch. (PR#19)
- return coordinates while training for debugging - returns status of the point such as revcomp,peak/nonpeak,coordinates etc

## [Unreleased] - 2022-01-15
- In `step6_train_chrombpnet_model` changed the file saving name bug  
- Changed scripts in `src/evaluatuion/make_bigwigs/` to be compatible with new dir structure

## [Unreleased] - 2022-01-13
- Added a note in README that users can continue to use the pre-trained bias model from chrombpnet-lite repo, if thats what they have been using till now
- Changed the default modisco crop setting in tutorial from 1000 to 500
- Fixed a bug in `src/evaluation/modisco/run_modisco.py` in directory creation
