# ChromBPNet: Deep learning models of base-resolution chromatin profiles
                This repo is under construction. Please check back in a week

- This repo contains code for the chrombpnet paper Anusri Pampari*, Anna Shcherbina*, Anshul Kundaje. (*authors contributed equally to the work)  
- Most have these scripts have been derived from an earlier repo (KerasAC) made by Anna Scherbina.
- General queries/thoughts have been addressed in the discussion section below.
- Please contact [Anusri Pampari] (\<first-name\>@stanford.edu)  for suggestions and comments. More instructions about reporting bugs detailed below.

## Quick Links

- [About](#chrombpnet)
- [Pipeline Overview](#pipeline-overview)
- [Training ChromBPNet](#training-chrombpent)
- [Testing ChromBpet ](#testing-chrombpnet)
- [Validation Experiments](#validation-experiemnts)
- [Helper Scripts](#helper-scripts)
- [Discussion](#discussion)
- [Report a bug](#chrombpnet-bugs)

##  ChromBPNet


##  Pipeline Overview

Refer to `testing.sh` for an example script to run the entire pipeline from end-to-end.


##  Training ChromBPNet


##  Testing ChromBPNet


##  Validation Experiments

- Interpretations
- Modisco/TomTom
- Marginal footprinting
- Invivo footprinting
- Variant Effect Prediction

##  Helper Scripts

- Preprocessing: Refer to the folder `src/utils/preprocessing/` for the following scrips (1) to convert bams to bigwigs (2) To generate tiledb database.
- Importance of correct shift: To understand why we are using +4/-4 shift for ATAC instead of the te popularly used +4/-5 shift refer to the scripts here - `src/utils/preprocessing/analysis` for the following scripts (1) Build PWM, (2) Train a simple model on random regions 
- Generate GC-Matched Negatives: To generate negatives that gc-match a given forground peak set refer to the folder `src/utils/preprocessing/make_gc_matched_negatives` 
- Choosing Hyperparameters: 
- Data generator examples: To get an understanding of the generators used in chrombpnet refer to the scripts here `src/utils/data_generator_examples`.

##  Discusssion

#### 

##  ChromBPNet bugs
