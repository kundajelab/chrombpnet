# Changelog
## [Unreleased]
- Updated the dev .yml file to only push when a push is made to the dev branch (after it got merged to main, it was initially pushing with changes to main as well)

## [0.2.1] - 2021-11-29
### Added
- Dockerfile to configure official Docker image on DockerHub
- Github workflow to automatically update Docker image whenever a new version is released (ie. whenever a tag named v* is pushed)
- Another Github workflow to automatically update a separate docker image (kundajelab/chrombpnet-lite:dev) whenever a push is made to dev

### Changed
- In scripts/predict_to_bigwig.py, removed "import keras" which was causing issues. 
- Updated requirements.txt to use Modisco 0.5.16 (the newest version) increased the max-seqlets in argmanager.py to 1 million


## [0.2.0] - 2021-10-19

Many changes. Upgraded to TF2 + changes to architecture.

### Added
- Mechanism to internally sample negatives provided a specific `--negative-sampling-ratio`, along with an option for `--no-negative-sampling`
- Support for models for which input length != output length (to accomodate valid padding arch)
- Specify threshold factor for bias model training filtering of non-peaks based on counts
- Abillity to specifiy profile and/or counts for interpretation (both by default)
- Debug chr option for interpretation
- Pretrained bias model for ATAC
- Notebooks to inspect trianed bias and ChromBPNet models

### Changed
- Changes for TF2, tested on v2.5.1, note v2.6.0 has memory leak that fails during chrombpnet training after a few epochs
- Model updated:
  - Now using valid padding for dilated convs, 2114->1000, 512 filters and 8 dilated convs default
- Command line arguments changed to specify `output_prefix` in many cases instead of `output_dir` for added flexibility
- Test chr for metrics has no default
- Lower memory requirements for training by switching to 8 bit representations for one-hot encoded seq
- Checks to not train if model with given prefix already exists

## [0.1.1]
### Added
- First version working with TF1.14, Keras 2.2.4.
- Training, metrics, interpretation and modsico supported and functional
- Architecture with same padding, same input and output sizes

--- 
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).


[unreleased]: https://github.com/kundajelab/chrombpnet-lite/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/kundajelab/chrombpnet-lite/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/kundajelab/chrombpnet-lite/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/kundajelab/chrombpnet-lite/releases/tag/v0.1.1
