source("00_init.R")
source("utils_plotting.R")

data_scaled <- readRDS("data/processed/metabolites_clean_no_outliers.rds")
k_range <- 2:16

k = 13 # change k if one wants to look at a different cluster distribution for HC when comparing linkage methods

# ==================================
# Hierarchical Clustering
# ==================================
## creates dendrograms to visually inspect clustering with HC and table with cluster assignments for a chosen k
compare_hc_linkage <- function (data,
                                k,
                                dendro_file = NULL,
                                cluster_file = NULL) {
  linkage_methods <- c("complete",
                       "ward.D2",
                       "single",
                       "average",
                       "mcquitty",
                       "centroid")
  hc_list <- list()
  cluster_tables <- data.frame(Cluster = 1:k)
  
  for (method in linkage_methods) {
    hc <- hclust(dist(data), method = method)
    hc_list[[method]] <- hc
    clusters <- cutree(hc, k = k)
    cluster_counts <- as.data.frame(table(clusters))
    colnames(cluster_counts) <- c("Cluster", method)
    
    cluster_tables <- merge(cluster_tables, cluster_counts, by = "Cluster")
  }
  
  cluster_tables <- cluster_tables[, -1]
  
  if (!is.null (dendro_file)) {
    pdf(dendro_file, width = 12, height = 6)
    par(mfrow = c(1, 2))
    for (method in linkage_methods) {
      plot(hc_list[[method]], main = paste0("Dendrogram (", method, ")"))
    }
    dev.off()
  }
  
  if (!is.null (cluster_file)) {
    write.csv(cluster_tables, cluster_file, row.names = TRUE)
  }
}

compare_hc_linkage(
  data_scaled,
  k,
  "output/plots/hc_linkage_dendrograms.pdf",
  "output/hc_linkage_cluster_tables.csv"
)

## decision for ward.D2

# ==================================
# k-means Clustering
# ==================================
### try conventional methods of finding k
# Elbow method
fviz_nbclust(
  data_scaled,
  kmeans,
  k.max = 16,
  method = "wss",
  print.summary = T
) +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(data_scaled, kmeans, k.max = 16, method = "silhouette") +
  labs(subtitle = "Silhouette method")

# # Gap statistic
set.seed(4242)
custom_kmeans <- function(data_scaled, centers, nstart) {
  kmeans(
    data_scaled,
    centers = centers,
    nstart = 100,
    iter.max = 100
  )
}

fviz_nbclust(
  data_scaled,
  custom_kmeans,
  nstart = 100,
  method = "gap_stat",
  nboot = 500
) +
  labs(subtitle = "Gap statistic method")

## Nb clust
NbClust(data_scaled, distance = "euclidean", method = "kmeans")

# ==================================
# Gaussian Mixture Models
# ==================================
# to try all models uncomment this section:
#fit <- Mclust(data_scaled, G = k_range)

# after trying all models, two were really close to each other and were looked at further to see differences
fit <- Mclust(data_scaled,
              G = k_range,
              modelNames = c("EII", "EEI"))

## to save BIC plot uncomment:
#pdf("new/output/plots/BIC_plot.pdf")

plot(fit, what = "BIC") ## ---> EEI 6, dann EEI 7, dann EEI 5
dev.off()
fit$BIC


# ==================================
# SOM (self organizing map) https://www.r-bloggers.com/2021/04/self-organizing-maps-in-r-supervised-vs-unsupervised/
# ==================================
data_matrix <- as.matrix(data_scaled)

### parameter tuning som
xdim_values <- c(5, 7, 10, 12)
ydim_values <- c(5, 7, 10, 12)
topo_values <- c("rectangular", "hexagonal")
radius_values <- c(10, 7, 5, 3, 2) # Beispielwerte für den Radius
rlen_values <- c(300, 500, 800, 1000)

# data must be a data.matrix
parameter_tuning_som <- function(data,
                                 xdim_values,
                                 ydim_values,
                                 topo_values,
                                 radius_values,
                                 rlen_values,
                                 save_pdf = TRUE) {
  results <- list()
  result_index <- 1
  
  # go through all parameter options
  for (xdim in xdim_values) {
    for (ydim in ydim_values) {
      for (topo in topo_values) {
        for (radius in radius_values) {
          for (rlen in rlen_values) {
            # create grid
            grid <- somgrid(xdim = xdim,
                            ydim = ydim,
                            topo = topo)
            
            # train SOM
            set.seed(4242)
            som_model <- som(
              data,
              grid = grid,
              alpha = c(0.05, 0.01),
              radius = radius,
              rlen = rlen
            )
            
            clusters <- som_model$unit.classif
            
            qual <- somQuality(som_model, data)
            te <- qual$err.topo
            qe <- qual$err.quant
            
            # save results
            results[[result_index]] <- list(
              xdim = xdim,
              ydim = ydim,
              topo = topo,
              radius = radius,
              rlen = rlen,
              som_model = som_model,
              clusters = clusters,
              qe = qe,
              te = te,
              silhouette = mean(silhouette(clusters, dist(data))[, 3])
            )
            result_index <- result_index + 1
          }
        }
      }
    }
  }
  
  # Sort results and print best QE and best TE
  results_df <- data.frame(
    xdim = sapply(results, function(x)
      x$xdim),
    ydim = sapply(results, function(x)
      x$ydim),
    topo = sapply(results, function(x)
      x$topo),
    radius = sapply(results, function(x)
      x$radius),
    rlen = sapply(results, function(x)
      x$rlen),
    qe = sapply(results, function(x)
      x$qe),
    te = sapply(results, function(x)
      x$te)
  )
  
  top_10_qe <- head(results_df[order(results_df$qe), ], 10)
  cat("Top 10 Quantization Errors:\n")
  print(top_10_qe)
  
  top_10_te <- head(results_df[order(results_df$te), ], 10)
  cat("\nTop 10 Topographic Errors:\n")
  print(top_10_te)
  
  saveRDS(results_df, "output/SOM_parameter_tuning.RDS")
  return (results_df)
}

results_tuning_som <- parameter_tuning_som(data_matrix,
                                           xdim_values,
                                           ydim_values,
                                           topo_values,
                                           radius_values,
                                           rlen_values)

## if tuned once it can be loaded here for scatter plot:
#results_tuning_som <- readRDS("output/SOM_parameter_tuning.RDS")

scatter_plot_som_tuning(data_matrix, results_tuning_som, pdf_file = "output/plots/Scatter_plot_SOM_parameter_selection.pdf")

## manual tuning and inspection:
grid <- somgrid(xdim = 12, ydim = 12, topo = "hex") # xdim und ydim größe des gitters, topo= sechseckig oder rechteckig ("rect")
pca_result <- prcomp(data_scaled)
data_pca <- pca_result$x[, 1:100]
data_pca <- as.matrix(data_pca)

set.seed(4242)
map <- som(
  data_pca,
  grid = grid,
  rlen = 300,
  # = epochs
  alpha = c(0.01, 0.005),
  # = lernrate (start und endwerte)
  radius = 10
)

# print QE and TE
qual <- somQuality(map, data_pca)
te <- qual$err.topo
print(te)
qe <- qual$err.quant
print(qe)

plot(map, type = "counts")
plot(map, type = "dist.neighbours")


# ==================================
# Subspace Clustering
# ==================================
#options(java.parameters = "-Xmx4g")
#library(subspace)
#library(rJava)

#xi = in wieviele Intervalle wird jeder metabolit unterteilt
#tau = "anzahl an datenpunkten damit Zelle als dicht gilt" 40% bei 0.4

## manual tunining of CLIQUE to get better understanding
system.time({
  clique_result <- CLIQUE(data_scaled, xi = 7, tau = 0.4)
})
clique_result
clique_result <- clique_result[order(unlist(lapply(
  lapply(clique_result, "[[", "subspace"), "sum"
)), decreasing = TRUE)]

clusters_clique <- clique_result$clusters.subspace

## parameter tuning with different values for Xi and Tau
xi_values <- c(20, 30, 40)
#,5,4,3)

tau_values <- c(0.1, 0.05, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

results_list <- list()
res_counter <- 0

for (xi in xi_values) {
  for (tau in tau_values) {
    cat("\n---------------------------------------\n")
    cat("start CLIQUE with xi =", xi, ", tau =", tau, "...\n")
    
    # shows what maximal dimension was in which clusters were found and saves results if clusters in more than 3D were found
    dims_per_cluster <- sapply(clique_res, function(cl)
      sum(cl$subspace))
    max_dimensions <- max(dims_per_cluster)
    
    cat(
      "Clusters found:",
      length(clique_res),
      "maximal dimensions in which clusters were found:",
      max_dimensions,
      "\n"
    )
    
    if (any(dims_per_cluster >= 3)) {
      cat("There was at least one cluster found in more than 3 dimensions! save results ... \n")
      
      res_counter <- res_counter + 1
      results_list[[res_counter]] <- list(
        xi = xi,
        tau = tau,
        clique_result = clique_res,
        dims_per_cluster = dims_per_cluster
      )
    } else {
      cat("no clusters found in more than 3 dimensions. results not saved.\n")
    }
  }
}


cat("\n\n--- Summary ---\n")
cat("In total",
    res_counter,
    "parameter combinations created at least 3D-Cluster.\n")

### find out which metabolites where the dimensions in which clusters > 4 were found

# loop through results
for (i in 1:length(results_list)) {
  result <- results_list[[i]]
  xi <- result$xi
  tau <- result$tau
  clique_result <- result$clique_result
  dims_per_cluster <- result$dims_per_cluster
  
  cat("\n---------------------------------------\n")
  cat("Results for xi =", xi, ", tau =", tau, ":\n")
  
  # loop through clusters
  for (j in 1:length(clique_result)) {
    cluster <- clique_result[[j]]
    dimensionen <- which(cluster$subspace) # Dimensionen (Metaboliten), die TRUE sind
    number_dimensionen <- length(dimensionen)
    
    # only consider clusters with more than 4 dimensions
    if (number_dimensionen >= 4) {
      samples <- cluster$objects
      sample_names <- rownames(data_scaled)[samples]
      
      cat("  Cluster", j, ":\n")
      cat("   Dimensions (Metaboliten):",
          paste(colnames(data_scaled)[dimensionen], collapse = ", "),
          "\n")
      cat("   Number of Samples:", length(samples), "\n")
      cat("   Sample names:", paste(rownames(data_scaled)[samples], collapse = ", "), "\n")
    }
  }
}

#saveRDS(results_list, file = "output/results_list_6_bis_9.rds")
results_list <- readRDS("output/results_list_6_bis_9.rds")

# ==================================
# Consensus Clustering -> not used in final Bachelor Thesis
# ==================================
# #### kmeans ####
# results_consensus_km <- ConsensusClusterPlus(as.matrix(t(data_scaled)),
#                                              maxK = 15,
#                                              reps = 1000,
#                                              pItem = 0.8,
#                                              pFeature = 1,
#                                              seed = 4242,
#                                              clusterAlg = "km",
#                                              distance = "euclidean",
#                                              plot = "png"
# )
#
# icl = calcICL(results_consensus,title="ICL",plot="png")
# icl[["clusterConsensus"]]
#
# saveRDS(results_consensus_km, "results/results_consensus_km.RDS")
# results_consensus_km <- readRDS("results_consensus_km.RDS")
# table(results_consensus_km[[10]]$consensusClass)
# rm(consensus_km)
# consensus_km <- data.frame(matrix(nrow = 428, ncol = 0))
# for (i in 2:15) {
#   # Spaltenname mit "cons" + i erstellen
#   column_name <- paste0("k_", i)
#
#   # Werte aus results_consensus_km[[i]]$consensusClass in die Spalte schreiben
#   consensus_km[[column_name]] <- results_consensus_km[[i]]$consensusClass
# }
# consensus_km[["k_16"]] <- results_consensus_km[[15]]$consensusClass