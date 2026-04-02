time0 <- Sys.time()
source("Input/config.R")
dir.create(outputroute, showWarnings = FALSE)
rmarkdown::render("Code/protein_level.Rmd", output_file = paste0("../", outputroute, project, " protein_results.html"))
Sys.time() - time0