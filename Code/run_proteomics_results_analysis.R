
source("Code/config.R")
rmarkdown::render("Code/protein_level.Rmd", output_file = paste0("../", outpuroute, project, " protein_results.html"))
