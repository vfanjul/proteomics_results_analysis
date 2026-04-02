

# Setup -------------------------------------------------------------------
source("Input/config.R")
source("Code/functions_general.R")
get_libraries(libraries, map_org(species)$go)
source("Code/functions_protein_level.R")
get_libraries(libraries, bioc_libraries)

## Experimental design
exp_design <- fread(designroute
)[, label := as.character(label)
][, group := factor(group, levels = unique(group))
][, sample := paste(group, replicate)
]

groups <- levels(exp_design$group)
samples <- exp_design$sample
sample_values <- samples
if (rel_to != "control_mean") sample_values <- paste0("∆ ", samples)

ctrl_group <- as.character(exp_design[control == 1, group][1])
case_groups <- setdiff(groups, ctrl_group)

group_means <- paste0("Mean z ", groups)
ctrl_mean <- paste0("Mean z ", ctrl_group)
case_params <- setdiff(group_means, ctrl_mean)
if (rel_to != "control_mean") case_params <- gsub("Mean", "∆ mean", case_params)

case_p_values <- gsub(".* z ", "P value ", case_params)
case_t_stats <- gsub(".* z ", "t statistic ", case_params)
case_change <- gsub(".*z ", "Change ", case_params)

case_dys <- gsub(".*z ", "dys ", case_params)
case_up <- gsub(".*z ", "up ", case_params)
case_down <- gsub(".*z ", "down ", case_params)

samples_by_group <- split(exp_design$sample, exp_design$group)
size_groups <- unlist(lapply(samples_by_group, length))



# Data preparation --------------------------------------------------------
## Protein data
prot_data <- fread(protroute)
setnames(prot_data, exp_design$label, samples)

prot_data <- get_prot_cols(prot_data, prot_col, species)

prot_data <- add_group_means(prot_data, samples_by_group, group_means, 
                             case_params, ctrl_mean, rel_to)


## Artifact Removal
prot_data <- flag_artifacts(prot_data, prot_col, rm_artifacts, exclude_pattern)
prot_data[Artifact == TRUE, .(Protein_name)][order(Protein_name)]

prot_data <- flag_inconsistent(prot_data, samples_by_group[[ctrl_group]], rm_inconsistent, quant_lim)
prot_data[Inconsistent == TRUE & Artifact != TRUE, .(Protein_name, `SD control`)][order(Protein_name)]


## Statistical significance
fit_cont <- fit_limma(prot_data[Artifact + Inconsistent == 0 & get(np_col) > 1], samples, exp_design$group)

prot_data[Artifact + Inconsistent == 0 & get(np_col) > 1, 
          c(gsub(paste0("-?", ctrl_group, "-?"), "", fit_cont$dt.cols),
            gsub("P value", "t statistic", gsub(paste0("-?", ctrl_group, "-?"), "", fit_cont$dt.cols))) :=
            c(as.list(as.data.frame(fit_cont$p.value)), as.list(as.data.frame(fit_cont$t)))]

zlim_opt <- optimize_zlim(prot_data[Artifact + Inconsistent == 0 & get(np_col) > 1],
                          case_params,
                          case_p_values,
                          zlim = zlims,
                          alpha = alphas)
facets_volcano <- arrange_facets(length(case_params))

prot_data <- flag_z_changes(prot_data, case_params, case_change, zlim_opt$best_zlim)
prot_data <- setorderv(prot_data, c(case_change, case_params), order = -1)

prot_data2 <- copy(prot_data)[Artifact + Inconsistent == 0 & get(np_col) > 1]


# Data exploration --------------------------------------------------------

## Normality
norm_data <- melt(prot_data[Artifact + Inconsistent == 0], id.vars = "Protein", measure.vars = samples, variable.name = "Sample", value.name = "Zq")
boxplot(Zq ~ Sample, norm_data, ylim = quantile(norm_data$Zq, c(0.01, 0.99)), 
        border = exp_design$color, col = NULL)
abline(h = c(-0.5, 0, 0.5), lty = c(3, 2, 3), col = "gray60")
any(abs(boxplot(Zq ~ Sample, norm_data, plot = FALSE)$stats[3, ]) > norm_lim)

plot_facets(prot_data[Artifact + Inconsistent == 0], samples, "", plot_qq)


## Variability
plot_dendrogram(prot_data2, samples)

par(mfrow = c(3, 1))
for (i in 2:4) plot_pca(prot_data2, samples, exp_design$color, groups, 1, i)
# for (i in 2:4) plot_pca(prot_data[Artifact + Inconsistent == 0], samples, exp_design$color, groups, 1, i)
# for (i in 2:4) plot_pca(prot_data2, samples, exp_design$color, groups, 1, i, TRUE)
par(mfrow = c(1, 1))

plot_legend(groups, unique(exp_design[, color]))


## Correlations
plot_corrpairs(prot_data2, samples)


## Protein changes
plot_facets(prot_data2,
            case_params,
            case_p_values,
            plot_volcano,
            zlim = zlim_opt$best_zlim,
            alpha = zlim_opt$best_alpha)



# Differential expression analysis ----------------------------------------

## Barplot
plot_bars(prot_data2, groups, case_params,
          unique(exp_design[, color]), zlim_opt$best_zlim, prot_col)


## Euler
prot_data2[, (case_dys) := lapply(case_params, function(x) abs(get(x)) >= zlim_opt$best_zlim)]
prot_data2[, (case_up) := lapply(case_params, function(x) get(x) >= zlim_opt$best_zlim)]
prot_data2[, (case_down) := lapply(case_params, function(x) -get(x) >= zlim_opt$best_zlim)]

if (length(groups) > 2) plot_euler(prot_data2,
                                   case_dys, case_groups,
                                   unique(exp_design[control == 0, color]))

if (length(groups) > 2) plot_euler(prot_data2,
                                   case_up, case_groups,
                                   unique(exp_design[control == 0, color]))

if (length(groups) > 2) plot_euler(prot_data2,
                                   case_down, case_groups,
                                   unique(exp_design[control == 0, color]))


## Heatmap
prot_data3 <- prot_data2[prot_data2[, Reduce(`|`, Map(`!=`, .SD, 0)) , .SDcols = case_change]]
prot_data3 <- prot_data3[prot_data3[, Reduce(`|`, Map(`<`, .SD, 0.05)) , .SDcols = case_p_values]]
prot_data3[, prot_lab := paste(Accession, Protein_name)]

if (sat_lim_95) {
  sat_lim <- round(quantile(abs(unlist(prot_data3[, sample_values, with = FALSE])), 0.95))
}

heat_width <- boxplot(strwidth(prot_data3[, prot_lab], units = "inches"), plot = FALSE)$stats[5] + length(samples)/10
heat_height <- (prot_data3[, .N] + 1)/10


pdf("Protein heatmap.pdf", heat_width, heat_height)
heatmap.2(as.matrix(prot_data3[, sample_values, with = FALSE]),
          Rowv = NA, Colv = NA, dendrogram = "none",
          ColSideColors = exp_design$color,
          labRow = prot_data3[, prot_lab], labCol = "",
          margins = c(0, 0),
          lmat = rbind(c(5, 4, 0), c(0, 1, 0), c(3, 2, 0)),
          # lhei = c(0.001, 0.1, 10),
          lhei = c(lcm(0.01*2.54), lcm(0.09*2.54), lcm(prot_data3[, .N]/10*2.54)),
          lwid = c(0.1, length(samples)/2,
                   boxplot(strwidth(prot_data3[, prot_lab], units = "inches"), plot = FALSE)$stats[5])/10,
          # lwid = c(0.01, lcm(length(samples)/2),
          #          lcm(boxplot(strwidth(prot_data3[, prot_lab], units = "inches"), plot = FALSE)$stats[5]/10)),
          col = colorpanel(sat_lim*100 - 1, "dodgerblue", "white", "red"),
          breaks = seq(-sat_lim, sat_lim, length.out = sat_lim*100),
          rowsep = (1:nrow(prot_data3))[!duplicated(prot_data3[, paste(case_change), with = FALSE])][-1] - 1,
          scale = "none",
          trace = "none",
          key = FALSE)

dev.off()

plot_legend(groups, unique(exp_design[, color]))
plot_prot_heatmap_key(sat_lim)



## Export data
data_cols <- unique(c(prot_col, np_col, samples, group_means, case_params, 
                      case_p_values, case_t_stats, case_change))
if (rm_artifacts == TRUE) data_cols <- c(data_cols, "Artifact")
if (rm_inconsistent == TRUE) data_cols <- c(data_cols, "Inconsistent")

dir.create(outputroute, showWarnings = FALSE)

fwrite(prot_data[, ..data_cols],
       paste0(outputroute, project, " protein_results.csv"), sep = ";")

write_excel(prot_data[, ..data_cols],
            paste0(outputroute, project, " protein_results.xlsx"), "protein_results",
            c(samples, case_params), case_p_values, case_change, np_col)




# Gene set expression analysis --------------------------------------------

## GO
gsea_go <- get_gsea(prot_data2, case_t_stats, 
                    species = species)
gsea_go_sig <- get_gsea_sig(gsea_go, filter_col = filter_col, sig_lim = sig_lim)
View(gsea_go_sig@result)

gsea_go_long <- get_gsea_long(gsea_go_sig, prot_data2, 
                              case_t_stats, prot_col)

go_height <- (length(unique(gsea_go_sig@result$Description)) + 20)/10

plot_gsea_ridges(gsea_go_long, y_max = cat_max, sat_lim = sat_lim)

for (i in case_params) plot_gsea_cnet(get_gsea_group(gsea_go, prot_data2, i, 
                                                     filter_col = filter_col, sig_lim = sig_lim), sat_lim = sat_lim)


# KEGG
kegg_db <- get_kegg_db(species = species, 
                       l1_blacklist = l1_blacklist,
                       l2_blacklist = l2_blacklist,
                       l3_blacklist = l3_blacklist,
                       l2_whitelist = l2_whitelist,
                       l3_whitelist = l3_whitelist)
# View(kegg_db$path2category[, .N, .(L1, L2, exclude)])
# kegg_db$path2category[is.na(exclude), .N, .(L1, L2)]

# View(kegg_db$path2category)
# dcast(kegg_db$path2category[, .N, .(L1, exclude = !is.na(exclude))], L1 ~ exclude, value.var = "N")

gsea_kegg <- get_gsea(prot_data2, case_t_stats, 
                      source_db = "KEGG", name_col = "Entrez", species = species, kegg_db = kegg_db)
gsea_kegg_sig <- get_gsea_sig(gsea_kegg, filter_col = filter_col, sig_lim = sig_lim)
View(gsea_kegg_sig@result)

gsea_kegg_long <- get_gsea_long(gsea_kegg_sig, prot_data2, 
                                case_t_stats, prot_col, id_col = "Entrez")

kegg_height <- (length(unique(gsea_kegg_sig@result$Description)) + 20)/10

plot_gsea_ridges(gsea_kegg_long, y_max = cat_max, sat_lim = sat_lim)

for (i in case_params) plot_gsea_cnet(get_gsea_group(gsea_kegg, prot_data2, i, name_col = "Entrez",
                                                     filter_col = filter_col, sig_lim = sig_lim), sat_lim = sat_lim)



fwrite(gsea_go_long, 
       paste0(outputroute, project, " gsea_go_results.csv"), sep = ";")

write_excel(gsea_go_long, 
            paste0(outputroute, project, " gsea_go_results.xlsx"), "gsea_go_results", 
            "NES", c("pvalue", "p.adjust"), "value", "setSize")

fwrite(gsea_kegg_long, 
       paste0(outputroute, project, " gsea_kegg_results.csv"), sep = ";")

write_excel(gsea_kegg_long, 
            paste0(outputroute, project, " gsea_kegg_results.xlsx"), "gsea_kegg_results", 
            "NES", c("pvalue", "p.adjust"), "value", "setSize")

