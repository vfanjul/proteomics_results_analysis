#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021

#### Setup ####
## Study
author <- "Víctor Fanjul, PhD"
project <- "Premature vs Normal Aging"
species <- "Mus musculus" # Mus musculus Homo sapiens Sus scrofa
databases <- ""
rel_to <- "control_mean" # control_mean sample_mean
fdr_quant <- 0.01 # False discovery rate for quantification


## Statistical Parameters
norm_lim <- 0.5 # Threshold for deviations from normality.
zlims <- seq(1, 3, 0.1)
alphas <- c(0.01, 0.05)


## Plotting
sat_lim <- 3 # Saturation limit for plots
sat_lim_95 <- TRUE # Saturation limit is percentile 95
dpi <- 300 # Figure resolution
figw <- 7 # Figure width
figh <- 7 # Figure height
figs <- 0.9 # Figure scale
cat_max <- 50 # Max number of categories to plot


## Artifact & Inconsistent Protein Filters
rm_artifacts <- TRUE # Remove artifact proteins
rm_inconsistent <- TRUE # Remove proteins with highly inconsistent controls
exclude_pattern <- c("krt", "keratin", # Skin keratins from handlers
                     "trypsin", "prss", # Trypsin from digestion
                     "albumin", "\\btransferrin\\b", "haptoglobin", "ahsg", "serpina", 
                     "apolipoprotein (a|e)", "taurus", "bovine", # Bovine from medium
                     "hba", "Hbb", "hemoglobin", "myoglobin", "complement c", # Other serum
                     "ntamin") # Labeled as contaminants
exclude_text <- "keratin, serum, and trypsin contaminants"
quant_lim <- 0.99


## GSEA
filter_col <- "p.adjust" # Filter in GSEA. Either "NES", "pvalue, or "p.adjust"
sig_lim <- "0.05" # Absolute lim value for filter col. 

l1_blacklist <- "organismal|disease|drug"
l2_blacklist <- "virus|prokaryote|maps"
l3_blacklist <- "plant|-.*bacter|- other$|worm|fly|yeast|meiosis"
l2_whitelist <- "circulatory|cardiovascular|environmental|aging"
l3_whitelist <- "null"
kegg_curation_text <- "Drug development, organismal systems, human diseases, overview maps, non-mammalian biology, and other-tissue-specific categories were excluded. Exceptions were categories related to environmental adaptation, cardiovascular system/disease, and aging. "


## Name of cols in datasets
prot_col <- "Protein"
np_col <- "Np"


#### File settings ####
designroute <- "Input/design.csv"
protroute <- "Input/protein_results.csv"
catroute <- "Input/category_results.csv"
outputroute <- "Output/"
export_csv <- FALSE # Export results as uncompressed CSV files

