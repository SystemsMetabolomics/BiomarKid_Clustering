source("00_init.R")
source("utils_plotting.R")

data_scaled <- readRDS("data/processed/metabolites_clean_no_outliers.rds") 

# ==================================
# CHOOSE SETTINGS
# ==================================
k_range <- 2:16
bootstrapping <- FALSE ## turn to true if bootstrapping should be run  (takes some time)
create_SOM_plot <- FALSE # turn to true to plot cluster assignments of k=13 for SOM
# ==================================
#
# ==================================

## performs 4 evaluation metrics and creates table with values
evaluation_df <- function(data, df, cluster_assig, method, i) {
  dist_matrix <- dist(data)
  
  # Davies-Bouldin-Index: Kompaktheit und separation, je kleiner der wert desto besser
  davies_bouldin_index <- index.DB(data, cluster_assig)
  df$davies_bouldin[i - min(k_range) + 1] <- davies_bouldin_index[1]
  #print(paste("Davies-Bouldin-Index für k =", i, ":", davies_bouldin_index[1]))
  
  # Silhouette-Index: misst, wie gut ein Objekt zu seinem Cluster passt, Werte nahe 1 sind gut
  silhouette_avg <- silhouette(cluster_assig, dist_matrix)[, 3]
  df$silhouette[i - min(k_range) + 1] <- mean(silhouette_avg)
  #print(paste("Silhouette-Index für k =", i, ":", mean(silhouette_avg)))
  
  # Dunn-Index: Verhältnis von minimalem Abstand zwischen Clustern zu maximalem Durchmesser innerhalb eines Clusters, höhere Werte sind besser
  dunn_index <- dunn(dist_matrix, cluster_assig)
  df$dunn[i - min(k_range) + 1] <- dunn_index
  #print(paste("Dunn-Index für k =", i, ":", dunn_index))
  
  # Calinski-Harabasz-Index: Verhältnis von Between-Cluster-Varianz zu Within-Cluster-Varianz, höhere Werte sind besser
  #calinski_harabasz_index <- cluster.stats(data, cluster_assig)$ch
  calinski_harabasz_index <- calinhara(data_scaled, cluster_assig) #(fpc index)
  df$calinski_harabasz[i - min(k_range) + 1] <- calinski_harabasz_index
  #print(paste("Calinski-Harabasz-Index für k =", i, ":", calinski_harabasz_index))
  
  return(df)
}

# ==================================
# Hierarchical Clustering
# ==================================
## pearson
cor_matrix <- cor(t(data_scaled), method = "pearson")
pearson_dist_matrix <- 1 - cor_matrix
hc_pearson <- hclust(as.dist(pearson_dist_matrix), method = "ward.D2")

## euclidean
hc_euclidean <- hclust(dist(data_scaled), method = "ward.D2")

distances <- list(pearson = hc_pearson, euclidean = hc_euclidean)
evaluation_hc_distance <- function(distances, data) {
  evaluation_hc_results <- list()
  cluster_hc_list <- list()
  
  for (distance_name in names(distances)) {
    tmp_evaluation <- data.frame(k = k_range)
    cluster <- data.frame(row.names = rownames(data)) # Zeilennamen von data_scaled übernehmen
    hc <- distances[[distance_name]]
    
    for (i in k_range) {
      cluster_assignments <- as.data.frame(cutree(hc, k = i))  # Adjust 'k' or 'h' based on significant clusters/heights of significant cluster pairs
      colnames(cluster_assignments) <- "cluster"
      cluster[, paste0("k_", i)] <- cluster_assignments # Spalte hinzufügen mit Clusterzuweisungen
      tmp_evaluation <- evaluation_df(data,
                                      tmp_evaluation,
                                      cluster_assignments$cluster,
                                      "hc",
                                      i)
    }
    
    evaluation_hc_results[[distance_name]] <- tmp_evaluation
    cluster_hc_list[[distance_name]] <- cluster
    rm(cluster_assignments,cluster,tmp_evaluation,distance_name,i)
  }
  result <- list(evaluation = evaluation_hc_results, clusters = cluster_hc_list)
  return(result)
}

hc_output <- evaluation_hc_distance(distances, data_scaled)

pdf("output/plots/evaluation_plots_hc_pearson_vs_ward_D.pdf", width = 8, height = 6)
evaluation_plots(hc_output$evaluation, names(distances), "distances")
dev.off()

## decision for euclidean based on evaluation plots
saveRDS(hc_output$evaluation$euclidean,"output/evaluation_hc_results.RDS")
saveRDS(hc_output$clusters$euclidean, "output/cluster_hc_list.RDS")

# -----------------------------
# robustness
# -----------------------------
if (bootstrapping){
res_loop_boot_hc <- list()

for (i in k_range) {
  res_boot <- clusterboot(
    dist(data_scaled),
    k = i,
    clustermethod = disthclustCBI, #hc with dist matrix as input
    B = 1000,
    bootmethod = "boot",
    multipleboot = TRUE, 
    bscompare = TRUE,
    count = TRUE, # counting am screen
    method = "ward.D2",
    seed = 4242
  )
  print(res_boot)
  res_loop_boot_hc[[i]] <- res_boot
}

saveRDS(res_loop_boot_hc, "output/results_bootstrapping_hc.RDS")
}

## k = 13 has best jaccard

# ==================================
# k-means Clustering
# ==================================
seeds <- c(42, 4242, 568, 123, 999, 2323)
evaluation_km_results <- list()
cluster_km_liste <- list()

for (seed in seeds) {
  tmp_evaluation <- data.frame(k = k_range)
  cluster_km <- data.frame(row.names = rownames(data_scaled))
  
  for (i in k_range) {
    set.seed(seed)
    results <- kmeans(
      data_scaled,
      centers = i,
      nstart = 150, # verschiedene nstarts manual getestet, siehe supplementary files
      iter.max = 20
    ) 
    cluster_km[, paste0("k_", i)] <- results$cluster
    
    tmp_evaluation <- evaluation_df(data_scaled, tmp_evaluation, results$cluster, "km", i)
  }
  evaluation_km_results[[paste0("seed", as.character(seed))]] <- tmp_evaluation
  cluster_km_liste[[paste0("seed", as.character(seed))]] <- cluster_km
  rm(tmp_evaluation, results)
}

#### PLOTS
pdf("output/plots/evaluation_plots_km_nstart150.pdf", width = 8, height = 6)
evaluation_plots(evaluation_km_results, seeds, "seeds")
dev.off()

### decision for seed4242 based on evaluation plots and save results
saveRDS(evaluation_km_results$seed4242, "output/evaluation_km_results_nstart150.RDS")
saveRDS(cluster_km_liste$seed4242,"output/cluster_km_list_nstart150.RDS")

# -----------------------------
# robustness
# -----------------------------
#### Bootstrapping
if (bootstrapping){
res_loop_boot_km <- list()

for (i in k_range) {
  res_boot <- clusterboot(
    data_scaled,
    runs = 150, 
    krange = i,
    clustermethod = kmeansCBI, # kmeans
    B = 1000,
    bootmethod = "boot",
    multipleboot = TRUE, 
    bscompare = TRUE,
    count = TRUE, # counting am screen
    seed = 4242
  )
  print(res_boot)
  res_loop_boot_km[[i]] <- res_boot
}

saveRDS(res_loop_boot_km, "output/results_bootstrapping_kmeans.RDS")
}
## k = 7 has best jaccard

# ==================================
# Gaussian Mixture Models
# ==================================
models = c("EEI", "EII")
#"VEI", error
#"VVI", error
#"VVV") error

evaluation_gmm_results <- list()
cluster_gmm_liste <- list()

for (model in models) {
  tmp_evaluation <- data.frame(k = k_range)
  cluster_gmm <- data.frame(row.names = rownames(data_scaled))
  
  for (i in k_range) {
    set.seed(4242)
    fit <- Mclust(data_scaled, G = i, modelNames = model)
    results <- data.frame(cluster = fit$classification)
    cluster_gmm[, paste0("k_", i)] <- fit$classification
    
    tmp_evaluation <- evaluation_df(data_scaled, tmp_evaluation, results$cluster, "gmm", i)
  }
  evaluation_gmm_results[[model]] <- tmp_evaluation
  cluster_gmm_liste[[model]] <- cluster_gmm
  rm(cluster_gmm, tmp_evaluation, model, i)
}

#### PLOTS
pdf("output/plots/evaluation_plots_gmm.pdf", width = 8, height = 6)
evaluation_plots(evaluation_gmm_results, models, "models")
dev.off()

### decision for model based on evaluation plots
saveRDS(evaluation_gmm_results$EEI,"output/evaluation_gmm_results.RDS")
saveRDS(cluster_gmm_liste$EEI, "output/cluster_gmm_list.RDS")

# ==================================
# SOM (self organizing map)
# ==================================
grid <- somgrid(xdim = 12, ydim = 12, topo = "hex") # xdim und ydim größe des gitters, topo= sechseckig oder rechteckig ("rect")

pca_result <- prcomp(data_scaled)
data_pca_50 <- pca_result$x[, 1:50]
data_pca_50 <- as.matrix(data_pca_50)
data_pca_100 <- pca_result$x[, 1:100]
data_pca_100 <- as.matrix(data_pca_100)

pca_list <- list(
  data_matrix = as.matrix(data_scaled),
  data_pca_50 = data_pca_50,
  data_pca_100 = data_pca_100
)
evaluation_som_results <- list()
cluster_som_list <- list()

for (data in names(pca_list)) {
  set.seed(4242)
  map <- som(
    pca_list[[data]],
    grid = grid,
    rlen = 300,
    alpha = c(0.01, 0.005), # = lernrate (start und endwerte)
    radius = 10
  )
  
  tmp_evaluation <- data.frame(k = k_range)
  cluster_som <- data.frame(row.names = rownames(data_scaled))
  
  for (i in k_range) {
    superclust_hclust <- hclust(dist(map$codes[[1]]), "ward.D2")
    superclasses_hclust <- cutree(superclust_hclust, i)
    cluster_assignments <- as.data.frame(superclasses_hclust[map$unit.classif])
    colnames(cluster_assignments) <- "cluster"
    cluster_som[, paste0("k_", i)] <- cluster_assignments # Spalte hinzufügen mit Clusterzuweisungen
    
    tmp_evaluation <- evaluation_df(data_scaled,
                                    tmp_evaluation,
                                    cluster_assignments$cluster,
                                    "som",
                                    i)
  }
  evaluation_som_results[[paste0("som", as.character(data))]] <- tmp_evaluation
  cluster_som_list[[paste0("som", as.character(data))]] <- cluster_som
  #rm(cluster_assignments, cluster, tmp_evaluation, method, i)
}

pdf("output/plots/evaluation_plots_som.pdf", width = 8, height = 6)
evaluation_plots(evaluation_som_results, c("no pca", "50 Compartments (82%)", "100 Compartments (92%)"), "matrices")
dev.off()

### decision for model based on evaluation plots
saveRDS(evaluation_som_results$somdata_pca_50, "output/evaluation_som_results.RDS")
saveRDS(cluster_som_list$somdata_pca_50, "output/cluster_som_list.RDS")


### aweSOM Plot for k=13 that shows cluster assignments of grid cells (clustered with HC)
if (create_SOM_plot){
  set.seed(4242)
  results_som <- som(
    data_pca_50,
    grid = grid,
    rlen = 300,
    alpha = c(0.01, 0.005),
    radius = 10
  )
  
  superclust_hclust <- hclust(dist(results_som$codes[[1]]), "ward.D2")
  superclasses_hclust <- cutree(superclust_hclust, 13)
  
  aweSOMplot(som = map, type = "Cloud", data = data_scaled, 
             #variables = colnames(data_scaled), 
             superclass = superclasses_hclust)
}

# -----------------------------
# robustness
# -----------------------------
#### Bootstrapping
if (bootstrapping){
res_loop_boot_som <- list()

for (i in k_range){
  res_boot <- clusterboot(dist(map$codes[[1]]), 
                          k = i,
                          clustermethod = disthclustCBI,
                          B = 1000,            
                          bootmethod = "boot",
                          multipleboot = TRUE, 
                          bscompare = TRUE,
                          count = TRUE,
                          method = "ward.D2",
                          seed = 4242)
  print(res_boot)
  res_loop_boot_som[[i]] <- res_boot
} 

res_loop_boot_som[[12]]

saveRDS(res_loop_boot_som, "output/results_bootstrapping_som.RDS")
}

