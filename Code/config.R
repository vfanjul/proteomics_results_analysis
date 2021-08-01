#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021


#### Setup ####

zlim <- 1.5
meanzlim <- 2
rm_artifacts <- TRUE

prot_field <- "Protein"
np_field <- "Np"


#### File paths ####

designroute <- "Input/design.csv"
protroute <- "Input/Datasets/protein_results.csv"
catroute <- "Input/Datasets/category_results.csv"
outpuroute <- "Output/"