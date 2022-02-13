#### Proteomics results analysis
## Created by: Victor Fanjul, Aug-2021

#### Load functions ####

libraries <- c("car", "eulerr", "GGally", "gplots", "ggplot2", "ggnewscale", "ggridges")

bioc_libraries <- c("limma", "clusterProfiler", "enrichplot")
# org.Ss.eg.db org.Mm.eg.db


#### New functions ####

#' Add Group Means
#' 
#' @description Adds columns with mean z to each group in a proteomics dataset.
#' 
#' @param dt Proteomics dataset.
#' @param samples_by_group List of group sample vectors.
#' @param group_mean_cols Vector of names for each group z mean column.
#' @param case_delta_cols Vector of names for each case group delta z mean column.
#' @param ctrl_mean_col Name of the control z mean group column.
#' @param rel_to If not "control_mean" (default), delta z means are calculated.
#' 
#' @return Data table with additional columns.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-12-05)

add_group_means <- function(dt, samples_by_group, group_mean_cols, 
                            case_delta_cols, ctrl_mean_col,
                            rel_to = "control_mean") {
  
  dt <- as.data.table(dt)
  dt[, (group_mean_cols) := lapply(samples_by_group,
                                   function(x) rowMeans(dt[, x, with = FALSE]))]
  
  if (rel_to != "control_mean") {
    dt[, (case_delta_cols) := lapply(setdiff(group_mean_cols, ctrl_mean_col), 
                                     function(x) get(x) - get(ctrl_mean_col))
    # ][, `∆ mean z` := rowMeans(.SD), .SDcols = case_delta_cols]
    ]
  }
  dt
}



#' Arrange Facet Display
#' 
#' @description Determines number of facet rows and columns in a multifaceted plot.
#' 
#' @param facets Number of facets.
#' 
#' @return Vector with number of rows and columns.
#' 
#' @author Victor Fanjul (2021-08-29)

arrange_facets <- function(facets) {
  plot_rows <- ifelse(facets == 3, 1, round(sqrt(facets), 0))
  plot_cols <- ifelse(facets == 3, 3, ceiling(facets/plot_rows))
  c(plot_rows, plot_cols)
}



#' Fit Limma with Empirical Bayes Statistics
#' 
#' @description Fits a lmFit linear model and computes empirical Bayes statistics
#' for differential expression.
#' 
#' @param dt Proteomics dataset.
#' @param samples Vector of sample z column names.
#' @param sample_groups Vector of sample groups (matching "samples").
#' 
#' @return eBayes object.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' * limma (Bioconductor)
#' 
#' @author Victor Fanjul (2021-12-05)

fit_limma <- function(dt, samples, sample_groups) {
  
  dt <- as.data.table(dt)[, .SD, .SDcols = samples]
  samples <- gsub(" ", "_", samples)
  names(dt) <- samples
  sample_groups <- factor(gsub(" ", "_", sample_groups), 
                          levels = unique(gsub(" ", "_", sample_groups)))
  groups <- levels(sample_groups)
  
  model_matrix <- model.matrix(~0 + sample_groups)
  rownames(model_matrix) <- samples
  colnames(model_matrix) <- groups
  
  contrasts <- apply(combn(groups, 2), 2, paste, collapse = "-")
  cont_matrix <- makeContrasts(contrasts = contrasts, levels = model_matrix)
  
  fit <- lmFit(dt, model_matrix)
  fit <- eBayes(contrasts.fit(fit, cont_matrix))
  fit$dt.cols <- paste0("P value ", gsub("_", " ", colnames(fit$p.value)))
  fit
}



#' Flag Artifacts in Proteomics Dataset
#' 
#' @description Adds a column that flags protein artifacts in a proteomics dataset.
#' 
#' @param dt Proteomics dataset.
#' @param prot_field Protein column name.
#' @param rm_artifacts Whether to remove artifacts (TRUE/FALSE).
#' @param exclude_pattern Vector of artifact name patterns.
#' @param artifact_field Name for artifact column. Default is "Artifact".
#' 
#' @return Data table with additional column.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-12-05)

flag_artifacts <- function(dt, prot_field, rm_artifacts, exclude_pattern, 
                           artifact_field = "Artifact") {
  dt <- as.data.table(dt)[, (artifact_field) := FALSE]
  if (rm_artifacts) {
    dt[grep(paste0(exclude_pattern, collapse = "|"), 
            get(prot_field), ignore.case = TRUE), (artifact_field) := TRUE]
  }
}



#' Flag Changes in Z
#' 
#' @description Adds columns to a proteomics dataset that flag whether there are 
#' changes (and its sign) in z for each sample group.
#' 
#' @param dt Proteomics dataset.
#' @param case_param_cols Vector of case group (delta) mean z column names.
#' @param case_change_cols Vector of new case group change column names.
#' @param zlim Threshold for determining change in z (positive number).
#' 
#' @return Data table with additional columns.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (date)
# 2021-12-05
flag_z_changes <- function(dt, case_param_cols, case_change_cols, zlim) {
  
  dt <- as.data.table(dt)
  dt[, (case_change_cols) := lapply(case_param_cols, 
                                    function(x) sign(get(x)) * as.numeric(abs(get(x)) >= zlim))]
}



#' Get Protein Fields
#' 
#' @description Extracts information from the protein name field in a proteomics 
#' dataset and includes it as new columns.
#' 
#' @param dt Proteomics dataset.
#' @param prot_field Protein column name.
#' 
#' @return Data table with additional columns.
#' 
#' @details Extracts informaiton for Accession, Entry_name and Protein_name.
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-10-19)

get_prot_fields <- function(dt, prot_field) {
  
  dt <- as.data.table(dt)
  dt <- dt[, aux := gsub(" PE=.*", "", gsub("^.*?\\|", "", get(prot_field)))
  ][, Accession := gsub("\\|.*", "", aux)
  ][, aux := gsub("^[[:alnum:]]*\\|", "", aux)
  ][, Entry_name := gsub(" .*", "", aux)
  ][, aux := gsub("^.*? ", "", aux)
  ][grep(" GN=", aux), Gene := gsub(".* GN=", "", aux)
  ][, Protein_name := gsub(" OS=.*", "", aux)
  ][, aux := NULL
  ]
  dt
}



#' Optimize Threshold for Proteomics Changes
#' 
#' @description Determines the optimal threshold for proteomics changes considering
#' z scores from WSSP and p values from limma.
#' 
#' @param dt Proteomics dataset.
#' @param z_fields Vector of z score column names.
#' @param p_fields Vector of p value column names.
#' @param zlim Vector of z thresholds to asses. Default is 1.5.
#' @param alpha Vector of alpha values to asses. Default is 0.05.
#' @param maximize Metric to maximize in the optimization. Default is "f1".
#' 
#' @return List containing metrics for each condition, mean metrics, the best z 
#' threshold, the best alpha and the value of the maximized metric.
#' 
#' @details Assesses accuracy, precision, recall and f1 for each z, p value and 
#' sample combination. Determines the optimal z and alpha thresholds when 
#' maximizing the selected metric (f1 by default).
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-08-28 29)

optimize_zlim <- function(dt, z_fields, p_fields, 
                          zlim = 1.5, 
                          alpha = 0.05,
                          maximize = "f1") {
  
  dt <- as.data.table(dt)
  metrics <- data.table(variable = 0, zlim = 0, alpha = 0, tn = 0, fn = 0, fp = 0, tp = 0)[-1]
  
  for (var in 1:length(z_fields)) for (i in zlim) for (j in alpha) {
    metrics <- rbind(metrics, 
                     unlist(list(z_fields[var], i, j, 
                                 as.list(dt[, .(table(z = abs(get(z_fields[var])) >= i, 
                                                      p = get(p_fields[var]) < j))
                                 ][order(z, p)][, N])), recursive = FALSE))
    
  }
  
  metrics[, accuracy := (tp + tn) / (fp + fn + tp + tn)
  ][, precision := tp / (tp + fp)
  ][, recall := tp / (tp + fn)
  ][, f1 := 2 * precision * recall / (precision + recall)
  ]
  
  mean_metrics <- metrics[, .(metric = mean(get(maximize), na.rm = TRUE)),
                          .(zlim, alpha)][order(-metric)]
  
  list(metrics = metrics[order(variable, -get(maximize))],
       mean_metrics = mean_metrics,
       best_zlim = mean_metrics[1, zlim],
       best_alpha = mean_metrics[1, alpha],
       max_metric = mean_metrics[1, metric])
}



#' Plot Barplot
#' 
#' @description Makes a barplot for group comparison.
#' 
#' @param dt Proteomics dataset.
#' @param groups Vector of groups.
#' @param values Vector of value column names.
#' @param colors Vector of group colors.
#' @param zlim Threshold for determining change in z (positive number).
#' @param prot_field Protein column name.
#' @param xlab Label for x axis. Default is "Differentially expressed proteins".
#' 
#' @return Barplot
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-10-18)

plot_bars <- function(dt, groups, values, colors, zlim, prot_field,
                      xlab = "Differentially expressed proteins") {
  
  dt <- as.data.table(dt)
  dtn <- dt[, .N]
  
  dt <- melt(dt, id.vars = prot_field, measure.vars = values
  )[abs(value) >= zlim, .(.N, rel = round(.N/dtn*100, 1)), variable]
  
  groups <- rev(groups)[1:dt[, .N]]
  colors <- rev(colors)[1:dt[, .N]]
  
  
  par(xpd = TRUE, mar = c(3,3,1,1), mgp = c(1.5,0.25,0), tck = - 0.01)
  bp <- barplot(rev(dt[, N]), horiz = TRUE, border = NA, 
                xlim = c(0, 2*max(dt[, N])), 
                col = colors,
                xlab = xlab,
                names = groups)
  text(rev(dt[, N]) + 10, bp, labels = paste(rev(dt[, rel]), "%"), adj = c(0, 0.5))
  par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA)
}



#' Plot Correlation Pairs
#' 
#' @description Makes a matrix of correlation plots for each pair combination.
#' 
#' @param dt Proteomics dataset.
#' @param samples Vector of sample z column names.
#' @param color Color of regression lines. Default is "red".
#' 
#' @return Correlation plot.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' * GGally
#' 
#' @author Victor Fanjul (2021-10-18)

plot_corrpairs <- function(dt, samples, 
                           color = "red") {
  dt <- as.data.table(dt)
  ggpairs(dt[, samples, with = FALSE],
          lower = list(continuous = wrap("smooth", color = color, size = 0.5)))
}




#' Plot Dendrogram
#' 
#' @description Makes a dendrogram for comparing samples.
#' 
#' @param dt Proteomics dataset.
#' @param samples Vector of sample z column names.
#' @param method Agglomeration clustering method. Default is "average".
#' 
#' @return Dendrogram plot.
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-10-17)

plot_dendrogram <- function(dt, samples, 
                            method = "average") {
  
  dt <- as.data.table(dt)
  par(mgp = c(1.5, 0.25, 0), tck = - 0.01)
  plot(hclust(dist(t(dt[, .SD, .SDcols = samples])), method = method), 
       labels = samples, main = "", xlab = "", sub = "", hang = -1)
  par(mgp = c(3, 1, 0), tck = NA)
}



#' Plot Euler Diagram
#' 
#' @description Makes an Euler diagram to compare variables in a dataset.
#' 
#' @param dt Proteomics dataset.
#' @param vars Vector of variable column names in dataset.
#' @param labels Vector of labels for each variable.
#' @param colors Vector of colors for each variable.
#' 
#' @return Euler Diagram
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' * eulerr
#' 
#' @author Victor Fanjul (2021-10-18)

plot_euler <- function(dt, vars, labels, colors) {
  
  dt <- as.data.table(dt)
  plot(euler(dt[, vars, with = FALSE], 
             shape = ifelse(length(vars) > 3, "ellipse", "circle")), 
       labels = labels,
       fills = adjustcolor(colors, 0.6),
       edges = list(col = "white"),
       quantities = list(type = c("counts")))
}



#' Make Faceted Plot
#' 
#' @description Makes a plot composed of multiple facets.
#' 
#' @param dt Proteomics dataset.
#' @param x Vector of variables to assess in x axis.
#' @param y Vector of variables to assess in y axis.
#' @param FUN Plot function. Default is NULL.
#' @param ... Additional parameters.
#' 
#' @return Multi-faceted plot.
#' 
#' @author Victor Fanjul (2021-08-29)

plot_facets <- function(dt, x, y, FUN = NULL, ...) {
  facets <- length(x)
  par(mfrow = arrange_facets(facets))
  for (i in 1:facets) FUN(dt, x[i], y[i], ...)
  par(mfrow = c(1, 1))
}



#' Plot Legend
#' 
#' @description Make a graphical legend.
#'
#' @param groups Vector of groups.
#' @param colors Vector of colors.
#' 
#' @return Legend plot
#' 
#' @author Victor Fanjul (2021-11-20)

plot_legend <- function(groups, colors) {
  
  par(xpd = TRUE, mar = c(0.1, 0.1, 0.1, 0.1))
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(x = "center", legend = groups, col = colors, pch = 19, bty = "n", ncol = length(groups))
  
  par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))
}



#' Plot Principal Components Analysis
#' 
#' @description Makes principal components analysis (PCA) and plots 2 components.
#' 
#' @param dt Proteomics dataset.
#' @param samples Vector of sample z column names.
#' @param colors Vector of colors.
#' @param groups Vector of groups.
#' @param comp_x Component displayed in x axis. Default is 1.
#' @param comp_y Component displayed in y axis. Default is 2.
#' 
#' @return PCA plot
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-10-17)

plot_pca <- function(dt, samples, colors, groups, 
                     comp_x = 1, 
                     comp_y = 2) {
  
  dt <- as.data.table(dt)
  protpca <- prcomp(t(dt[, .SD, .SDcols = samples]), scale. = TRUE)
  
  par(xpd = TRUE, mar = c(3,3,1,1), mgp = c(1.5,0.25,0), tck = - 0.01)
  plot(protpca$x[, c(comp_x, comp_y)], col = colors, pch = 19)
  legend(x = "top", legend = groups, col = unique(colors), pch = 19, cex = 0.5, 
         bty = "n", ncol = length(groups), inset = c(0,-0.1))
  
  par(xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA)
}



#' Plot Proteomics Heatmap
#' 
#' @description Makes a heatmap of proteomics changes.
#' 
#' @param dt Proteomics dataset.
#' @param samples Vector of sample z column names.
#' @param change_cols Vector of group change column names.
#' @param sample_colors Vector of colors for each sample
#' @param label_col Row label column name.
#' @param sat_lim Color saturation z threshold. Default is 3.
#' @param change_colors Vector with color scale for z and change columns. 
#' Default is c("dodgerblue", "white", "red").
#' 
#' @return Heatmap
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' * gplots
#' 
#' @author Victor Fanjul (2021-12-06)

plot_prot_heatmap <- function(dt, samples, change_cols, sample_colors, label_col, 
                              sat_lim = 3, 
                              change_colors = c("dodgerblue", "white", "red")) {
  
  dt <- as.data.table(dt)
  
  heatmap.2(as.matrix(dt[, samples, with = FALSE]),
            Rowv = NA, Colv = NA, dendrogram = "none",
            ColSideColors = sample_colors,
            labRow = dt[, get(label_col)], labCol = "",
            margins = c(0, 0),
            lmat = rbind(c(5, 4, 0), c(0, 1, 0), c(3, 2, 0)), 
            lhei = c(lcm(0.01*2.54), lcm(0.09*2.54), lcm(dt[, .N]/10*2.54)), 
            lwid = c(0.1, length(samples)/2,
                     boxplot(strwidth(dt[, get(label_col)], units = "inches"), plot = FALSE)$stats[5])/10,
            col = colorpanel(sat_lim*100 - 1, change_colors[1], change_colors[2], change_colors[3]), 
            breaks = seq(-sat_lim, sat_lim, length.out = sat_lim*100), 
            rowsep = (1:nrow(dt))[!duplicated(dt[, paste(change_cols), with = FALSE])][-1] - 1,
            scale = "none",
            trace = "none",
            key = FALSE)
}



#' Plot Heatmap Key
#' 
#' @description Makes a plot with a heatmap color key
#' 
#' @param sat_lim Color saturation z threshold. Default is 3.
#' @param change_colors Vector with color scale for z and change columns. 
#' Default is c("dodgerblue", "white", "red").
#' 
#' @return Heatmap key plot
#' 
#' @details 
#' 
#' # Required libraries:
#' * gplots
#' 
#' @author Victor Fanjul (2021-11-20)

plot_prot_heatmap_key <- function(sat_lim = 3, 
                                  change_colors = c("dodgerblue", "white", "red")) {
  
  par(mar = c(1.8, 0.2, 0.2, 0.2), mgp = c(0.8, 0.5, 0), tck = - 0.3, cex.axis = 0.8)
  plot(seq(-sat_lim, sat_lim, length.out = sat_lim*100), rep(0, sat_lim*100), 
       type = "n", axes = FALSE, xlab = "", ylab = "")
  
  axis(side = 1, at = seq(-sat_lim, sat_lim), line = 0.2)
  abline(v = seq(-sat_lim, sat_lim, length.out = sat_lim*100), lwd = 3, 
         col = colorpanel(sat_lim*100, change_colors[1], change_colors[2], change_colors[3]))
  
  par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), tck = NA, cex.axis = 1)
}




#' Plot Quantile-Quantile
#' 
#' @description Makes a Q-Q plot
#' 
#' @param dt Proteomics dataset.
#' @param x Column name with variable to assess.
#' @param ... Additional parameters.
#' 
#' @return Q-Q plot
#' 
#' @details 
#' 
#' # Required libraries:
#' * data.table
#' * car
#' 
#' @author Victor Fanjul (2021-10-17)

plot_qq <- function(dt, x, ...) {
  
  dt <- as.data.table(dt)
  
  par(mar = c(4, 4, 2, 2))
  qqPlot(dt[, get(x)], 
         xlab = "Theoretical normal quantiles",
         ylab = paste(x, " quantiles"), 
         col.lines = "red", pch = 16, cex = 0.5, id = FALSE)
  
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}



#' Make Volcano Plot
#' 
#' @description Makes a volcano plot.
#' 
#' @param dt Proteomics dataset.
#' @param z_field Z score column name.
#' @param p_field P value column name.
#' @param zlim Z threshold Default is 1.5.
#' @param alpha Alpha threshold. Default is 0.05.
#' @param cex Point size. Default is 0.5.
#' @param xlim Vector of x axis limits. Default is c(-5, 5).
#' @param ylim Vector of y axis limits. Default is c(0, 5).
#' @param color Vector of point colors. Default is c("red", "dodgerblue", "limegreen", "grey40").
#' @param transp Vector of point transparencies. Default is c(0.8, 0.6, 0.6, 0.4).
#' 
#' @return Volcano plot.
#' 
#' @details Point color and transparency may be personalized in each quadrant.
#' Order in vector is as follows:
#' * True positives: significance given by both z and p values. Default is red.
#' * False positives: significance given by z. Default is dodgerblue.
#' * False negatives: significance given by p values. Default is limegreen.
#' * True negatives: non-significant. Default is grey40.
#' 
#' # Required libraries:
#' * data.table
#' 
#' @author Victor Fanjul (2021-08-28)

plot_volcano <- function(dt, z_field, p_field, 
                         zlim = 1.5, 
                         alpha = 0.05,
                         cex = 0.5,
                         xlim = c(-5, 5),
                         ylim = c(0, 5),
                         color = c("red", "dodgerblue", "limegreen", "grey40"),
                         transp = c(0.8, 0.6, 0.6, 0.4)) {
  
  dt <- as.data.table(dt)
  
  par(pch = 16, cex = 1)
  plot(-log10(get(p_field)) ~ get(z_field), dt,
       xlab = z_field, ylab = bquote(~-Log[10]~ "p value"), 
       xlim = xlim, ylim = ylim, col = 0, pch = 16)
  
  filter <- list(dt[abs(get(z_field)) >= zlim & get(p_field) < alpha],
                 dt[abs(get(z_field)) >= zlim & get(p_field) >= alpha],
                 dt[abs(get(z_field)) < zlim & get(p_field) < alpha],
                 dt[abs(get(z_field)) < zlim & get(p_field) >= alpha])
  
  for (i in 1:4)  points(-log10(get(p_field)) ~ get(z_field), 
                         filter[[i]], cex = cex, 
                         col = adjustcolor(color[i], alpha.f = transp[i]))
  
  abline(v = c(-zlim, zlim), col = "grey60")
  abline(h = -log10(alpha), col = "grey60")
}



