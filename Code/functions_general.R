#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021

#### Load functions ####

libraries <- c("data.table", 
               "openxlsx",
               "BiocManager")

bioc_libraries <- c("org.Mm.eg.db", 
                    "GOstats")


#### New functions ####

#' Get Libraries
#' 
#' @description Installs and loads all required libraries.
#' 
#' @param general_libraries Vector of libraries in CRAN.
#' @param bioc_libraries Vector of libraries in Bioconductor. Default is NULL.
#' 
#' @details BiocManager is installed and loaded if required before the installation
#' of Bioconductor libraries.
#' 
#' @author Victor Fanjul (2021-08-01)

get_libraries <- function(general_libraries, 
                          bioc_libraries = NULL) {
  
  if (!is.null(bioc_libraries)) {
    general_libraries <- unique(c(general_libraries, "BiocManager"))
  }
  
  get_repo_libraries(general_libraries)
  if (!is.null(bioc_libraries)) get_repo_libraries(bioc_libraries, TRUE)
}



#' Get Libraries from Repository
#' 
#' @description Installs and loads required libraries from specific repository 
#' (CRAN or Bioconductor).
#' 
#' @param libraries Vector of libraries.
#' @param bioconductor Whether libraries are from Bioconductor. Default is FALSE.
#' 
#' @author Victor Fanjul (2021-08-01)

get_repo_libraries <- function(libraries, 
                               bioconductor = FALSE) {
  
  new_lbs <- setdiff(libraries, rownames(installed.packages()))
  unloaded_lbs <- setdiff(libraries, (.packages()))
  
  if (length(new_lbs) > 0) if (bioconductor) {
    BiocManager::install(new_lbs)
  } else install.packages(new_lbs)
  
  if (length(unloaded_lbs) > 0) lapply(unloaded_lbs, 
                                       library, 
                                       character.only = TRUE)
  
}



#' Write Excel File
#' 
#' @description Exports proteomics dataset to Excel with conditional formatting.
#' 
#' @param dt Proteomics dataset.
#' @param file File path and name.
#' @param worksheet Name of Excel sheet.
#' @param sample_cols Vector of z column names.
#' @param p_value_cols Vector of p value column names.
#' @param change_cols Vector of change column names.
#' @param np_col Number of peptides (Np) column.
#' @param sat_lim Color saturation z threshold. Default is 3.
#' @param change_colors Vector with color scale for z and change columns. 
#' Default is c("dodgerblue", "white", "red").
#' @param enhance_colors Vector with color scale for Np and p value columns.
#'  Default is c("limegreen", "white").
#' 
#' @details 
#' 
#' # Required libraries:
#' * openxlsx
#' 
#' @author Victor Fanjul (2021-12-05)

write_excel <- function(dt, file, worksheet, sample_cols, p_value_cols, 
                        change_cols, np_col, 
                        sat_lim = 3,
                        change_colors = c("dodgerblue", "white", "red"),
                        enhance_colors = c("limegreen", "white")) {
  
  data_fields <- names(dt)
  wb <- createWorkbook()
  addWorksheet(wb, worksheet)
  writeData(wb, worksheet, dt, withFilter = TRUE)
  
  conditionalFormatting(wb, worksheet,
                        cols = which(data_fields %in% sample_cols), 
                        rows = 1:nrow(dt) + 1,
                        style = change_colors,
                        rule = c(-sat_lim, 0, sat_lim),
                        type = "colourScale")
  
  conditionalFormatting(wb, worksheet,
                        cols = which(data_fields %in% p_value_cols), 
                        rows = 1:nrow(dt) + 1,
                        style = enhance_colors,
                        rule = c(0, 0.05),
                        type = "colourScale")
  
  conditionalFormatting(wb, worksheet,
                        cols = which(data_fields %in% change_cols), 
                        rows = 1:nrow(dt) + 1,
                        style = change_colors,
                        rule = c(-1, 0, 1),
                        type = "colourScale")
  
  conditionalFormatting(wb, worksheet,
                        cols = which(data_fields %in% np_col), 
                        rows = 1:nrow(dt) + 1,
                        style = enhance_colors,
                        rule = c(1, 2),
                        type = "colourScale")
  
  saveWorkbook(wb, file, TRUE)
}



