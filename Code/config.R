#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021


#### Setup ####

## Study
author <- "Víctor Fanjul, PhD"
project <- "Premature vs Normal Aging"
species <- "Mus musculus" # Mus musculus Homo sapiens Sus scrofa
databases <- ""
rel_to <- "sample_mean" # control_mean sample_mean
fdr_quant <- 0.01 # False discovery rate for quantification

## Parameters
zlims <- seq(1, 3, 0.1)
alphas <- c(0.01, 0.05)
sat_lim <- 3 # Saturation limit for plots

## Artifact Proteins
rm_artifacts <- TRUE # Remove artifact proteins
exclude_pattern <- c("krt", "keratin", "hba", "Hbb", "hemoglobin", "myoglobin", "haptoglobin", "albumin", "trypsin", "serpin", "complement c")
exclude_text <- "keratin contaminants, trypsin, and several major serum proteins (albumin, globins, serpins and complement factors)"

## Name of fields in datasets
prot_field <- "Protein"
np_field <- "Np"


#### File paths ####

designroute <- "Input/design.csv"
protroute <- "Input/protein_results.csv"
catroute <- "Input/category_results.csv"
outpuroute <- "Output/"
