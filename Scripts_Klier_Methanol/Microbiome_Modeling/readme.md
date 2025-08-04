### Scripts for the modeling part of "Integrating population-based metabolomics with computational microbiome modelling identifies methanol as a urinary biomarker for protective diet-microbiome-host interactions"


This repository contains all the necessary files for microbiome community modeling of the manuscript.

## Files

### `normalizedCoverage.csv`
Contains abundance table mapped onto AGORA2

### `run_mgPipe.m`
Starts the Microbiome Modeling Toolbox

### `run_MicrobePredictions.m`
Runs the predictMicrobeContributions function after creating microbiome community models with run_mgPipe.m

### `ReadOutMethanolStrains.m`
Identifies which strains in AGORA2 can secrete methanol

## Folders

### `/Diet_constraints`
Contains the diet constraints for the used diet constraints of the Average European Diet and additionally files of dient constraints
with added pectin/xylan constraints

