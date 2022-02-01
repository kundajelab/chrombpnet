# Changelog

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
