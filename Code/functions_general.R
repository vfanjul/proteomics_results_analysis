#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021

#### Load functions ####

libraries <- c("data.table", 
               "BiocManager")

bioc_libraries <- c("org.Mm.eg.db", 
                    "GOstats")

get_repo_libraries <- function(libraries, bioconductor = FALSE) {
  
  new_lbs <- setdiff(libraries, rownames(installed.packages()))
  unloaded_lbs <- setdiff(libraries, (.packages()))
  
  if (length(new_lbs) > 0) if (bioconductor) {
    BiocManager::install(new_lbs)
  } else install.packages(new_lbs)
  
  if (length(unloaded_lbs) > 0) lapply(unloaded_lbs, 
                                       library, 
                                       character.only = TRUE)
  
}


get_libraries <- function(general_libraries, 
                          bioc_libraries = NULL,
                          print_loaded = FALSE
                          ) {
  
  if (!is.null(bioc_libraries)) {
    general_libraries <- unique(c(general_libraries, "BiocManager"))
  }
  
  get_repo_libraries(general_libraries)
  if (!is.null(bioc_libraries)) get_repo_libraries(bioc_libraries, TRUE)
  
  if(print_loaded) {
    message("Loaded libraries are:")
    print((.packages()))
  }
}

get_libraries(libraries, bioc_libraries)



#### New functions ####