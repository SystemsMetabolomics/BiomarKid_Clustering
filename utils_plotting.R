source("utils_metabolite_relevance.R")


# ==================================
# See Data distribution (skewed or not) per metabolite and save as pdf
# ==================================
plot_histograms <- function(data, pdf_file = NULL, bins = 30, fill_color = "grey", border_color = "black") {
  if (!is.null(pdf_file)) {
    pdf(pdf_file)
  }
  
  # loop through metabolites and create histogramm
  for (metabolite in colnames(data)[-1]) {
    p_hist <- ggplot(data, aes_string(x = metabolite)) +
      geom_histogram(bins = bins, fill = fill_color, color = border_color, alpha = 0.7) +
      ggtitle(paste("Histogram for", metabolite)) +
      theme_minimal()
    
    print(p_hist)
    
  }
  if (!is.null(pdf_file)) {
    dev.off()
  }
}

# ==================================
# Plot PCA with marked outliers (outliers previously calculated with mahalanobis distance)
# ==================================
plot_pca_outliers <- function(pca_scores, outliers, pdf_file = NULL){
  if (!is.null(pdf_file)) {
    pdf(pdf_file)
  }
  
  # plotting of first 2 compartments of PCA
  plot(pca_scores[, 1], pca_scores[,2], main = "PCA with marked outliers", xlab = "PC1", ylab = "PC2", pch = 16) # pc = plotting character: 16 = filled dot
  
  # mark and label outliers red
  points(pca_scores[outliers, 1], pca_scores[outliers, 2], col = "red", pch = 19) # pch = 19: larger filled dot than 16
  text(pca_scores[outliers, 1], pca_scores[outliers, 2], labels = rownames(pca_scores)[outliers], col = "red", pos = 3) # pos = 3: above dot

if (!is.null(pdf_file)) {
  dev.off()
}
}

# ==================================
# Scatter Plot SOM tuning
# ==================================
scatter_plot_som_tuning <- function(data, results_tuning, pdf_file = NULL){
  if(!is.null(pdf_file)){
    pdf(pdf_file, width = 8, height = 6)
    par(mfrow = c(1, 1))
    # create scatterplot
    plot(results_tuning$qe, results_tuning$te,
         xlab = "Quantization Error (QE)",
         ylab = "Topographic Error (TE)",
         main = "QE vs. TE for SOM Parameter Selection")
    
    # labels of data points
    text(results_tuning$qe, results_tuning$te,
         labels = paste(results_tuning$xdim,
                        results_tuning$ydim,
                        results_tuning$topo,
                        results_tuning$radius,
                        results_tuning$rlen),
         pos = 4, cex = 0.5) # pos = 4 platziert die Beschriftungen rechts vom Punkt
    dev.off()
  }
}

# ==================================
# Internal Evaluation Plots
# ==================================
evaluation_plots <- function(evaluation_results, names_comparisons, comparison = as.character()){
  par(mfrow = c(2, 2))
  evaluation <- colnames(evaluation_results[[1]])[-1] #colnames of evaluation_df without k (only evaluation metrics)
  
  for (metric in evaluation){
    # calculate range of values for plot dimensions
    values <- numeric()
    for (df in evaluation_results){
      values <- c(values, df[[metric]]) 
    }
    value_range <- range(values)
    
    # plot first line
    plot(k_range, #k_values
         evaluation_results[[1]][[metric]], # each metric in the first dataframe to start first line with
         type = "l", # line plot
         main = paste(metric, "for different", comparison), # title of plots
         xlab = "k",
         ylab = metric,
         ylim = value_range # as calculated above
    )
    
    # plot other lines (for dataframes 2:last)
    for (i in 2:length(evaluation_results)){
      lines(k_range,
            evaluation_results[[i]][[metric]],
            col = i # for each new dataframe a new color
      )
      
      legend("topright", legend = names_comparisons, col = 1:length(names_comparisons), lty = 1, cex = 0.5)  # lty = 1 : whole line
    }
  }
}


# ==================================
# ari plots different methods (overlap clusters across methods)
# ==================================
ari_plots_methods <- function(methods){
  # num_k shows how many k's were "tested" in the methods. they should have same k (columns) because here just first one is taken as reference
  num_k <- ncol(methods[[1]])
  ari_results <- data.frame()
  
  # Comparison of all methodpairs for each k
  ## first loop: through each k
  for (k in 1:num_k) {
    ## second loop: through each method (j = first method of pair to be compared)
    for (j in 1:length(names(methods))) {
      # if statement, so last object is not compared to itself
      if (j + 1 <= length(names(methods))) {
        ## third loop: through each method (l = second method of pair to be compared)
        for (l in (j + 1):length(names(methods))) {
          # if statement to make sure objects are not compared to themselves
          if (j != l) {
            ari <- adjustedRandIndex(methods[[j]][[k]], methods[[l]][[k]])
            ari_results <- rbind(
              ari_results,
              data.frame(
                k = k + 1, # because k starts at 1 but cluster start at 2
                method1 = names(methods)[j],
                method2 = names(methods)[l],
                ARI = ari
              )
            )
          }
        }
      }
    }
  }
  
  # labels for plot
  ari_results$plotlabels <- paste(ari_results$method1, "vs.", ari_results$method2)
  
  # create plot with ggplot
  ggplot(ari_results,
         aes(
           x = k,
           y = ARI,
           color = plotlabels,
           group = plotlabels
         )) +
    geom_line(size = 0.5) +   # lines
    geom_point(size = 1) +  # dot for each combination
    labs(
      title = "Cluster Overlap across Methods with ARI",
      x = "Number of Cluster (k)",
      y = "ARI",
      color = "methods"
    ) 
}


# ==================================
# alluvial for different methods
# ==================================
alluvial_plot <- function (cluster_hc, cluster_km, cluster_gmm, cluster_som, k){
  ## create frequency table for given k 
  freq_table <- as.data.frame(table(
    cluster_hc[[paste0("k_", k)]],
    cluster_km[[paste0("k_", k)]],
    cluster_som[[paste0("k_", k)]],
    cluster_gmm[[paste0("k_", k)]]
  ))
  
  # add colnames
  colnames(freq_table) <- c("hc", "km", "som","gmm", "Freq")
  # only consider combinations with at least one object in them
  freq_table <- subset(freq_table, Freq > 0)
  
  p <- ggplot(data = freq_table,
              aes(axis1 = hc,  
                  axis2 = km, 
                  axis3 = som,
                  axis4 = gmm,
                  y = Freq)) +
    geom_alluvium(aes(fill = factor(km)), width = 1/2) +
    geom_stratum(width = 1/5) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("hc", "kmeans", "som", "gmm"), expand = c(0.15, 0.05)) +
    ggtitle(paste0("Alluvial Plot for different clustering methods (k = ",k , ")")) 
  #+ theme_minimal()
  
  return(p)
}


# ==================================
# Heatmaps
# ==================================
heatmap_every_sample_by_rows <- function (data, cluster_assignments, clip = 3, method, k){
  # clip data for heatmaply
  clipped_data <- data
  clipped_data[clipped_data > clip] <- clip
  clipped_data[clipped_data < -clip] <- -clip
  
  sorted_cluster <- sort_cluster(clipped_data, cluster_assignments, k)
  
  cluster_colors <- brewer.pal(k,"Set3")
  row_colors <- cluster_colors[as.factor(sorted_cluster$cluster)]
  sorted_cluster <- sorted_cluster[, !(names(sorted_cluster) == "cluster")] ## column "cluster" wieder entfernen bevor heatmap erstellt wird
  
  col<- colorRampPalette(c("red", "white", "blue"))(100)
  
  p <- heatmaply(
    sorted_cluster, 
    Rowv = FALSE,   # ordnet reihenfolge nicht selbst an (hab ich vorher ja schon angeordnet)
    xlab = "Metabolites",
    ylab = "Samples", 
    main = paste("Interactive Heatmap", method),
    col = col,
    RowSideColors = row_colors,
    showticklabels = c(TRUE, FALSE),
    fontsize_col = 4,
    file = paste0("output/plots/heatmaply_", method, "new_", k, ".html")
  )
  # Gesamte Legende ausblenden:
  p <- layout(p, showlegend = FALSE)
  p
}

plot_cluster_circos <- function(cluster_means, row_order, breaks = c(-3,0,3)){
  t_cluster_means <- t(scale(cluster_means)) 
  col_fun <- colorRamp2(c(-3,0,3),c("blue", "white", "red"))
  
  # initialize circos
  circos.clear()
  circos.par(start.degree = 60, gap.degree = 30)
  
  circos.heatmap(t_cluster_means[row_order, , drop = FALSE],
                 col = col_fun,
                 cluster = FALSE,
                 rownames.side = "outside",
                 rownames.cex = 0.35,
                 track.height = 0.6
  )
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 1) { # the last sector
      cn = rev(rownames(cluster_means))
      n = length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                  1:n - 0.5, cn, 
                  cex = 0.5, adj = c(0, 0.5), facing = "inside")
    }
  }, bg.border = NA
  )
  circos.clear()
}
