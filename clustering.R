source("data_preprocessing.R")
library(readxl)
library(writexl)
library(dplyr)
library(magrittr)
library(tidyr)
library(patchwork) # needed to combine plots
library(factoextra) # for PCA (screeplot)
library(ggplot2)
library(openxlsx)

data_scaled <- readRDS("data_scaled_and_adjusted.rds")

# ==================================
# SGI (Hierarchical Clustering) https://github.com/krumsieklab/sgi/tree/master
# ==================================
require(devtools)
#devtools::install_github(repo="krumsieklab/sgi", subdir="sgi")
library(sgi)

## format data
rownames(selected_phenotypes) <- selected_phenotypes$Sample
selected_phenotypes <- selected_phenotypes[,-1]

##### pearson correlation #####
#cor_test <- as.data.frame(data_log)
#rownames(cor_test) <- cor_test$Sample
#cor_test <- cor_test[,-1]
#cor_matrix <- cor(t(cor_test))  # Transpose data if clustering rows

# Convert correlation to distance
#cor_dist <- as.dist(1 - cor_matrix)

##### #####

## different hc methods (when done on raw data one can see clear outliers when comparing the plots)
hc_complete = hclust(dist(data_scaled), method = "complete")
hc_ward2 = hclust(dist(data_scaled), method = "ward.D2")
hc_single = hclust(dist(data_scaled), method = "single")
#hc_pear <- hclust(cor_dist, method = "ward.D2")

## save 3 plots
#png("combined_plots_rem_outliers3_wind_cadj.png", width = 1200, height = 600)
#par(mfrow = c(1, 3))  # 1 row, 3 columns
plot(hc_ward2, main = "Ward's Method", sub = "", xlab = "")
#plot(hc_complete, main = "Complete Method", sub = "", xlab = "")
#plot(hc_single, main = "Single Method", sub = "", xlab = "")
#plot(hc_pear, main = "Pearson Correlation Distance", sub = "", xlab = "")
#dev.off()

## chose method
hc <- hc_ward2

if (!isTRUE(all.equal(hc$labels, rownames(selected_phenotypes)))) {
  stop("Identifiers in hc$labels and rownames(selected_phenotypes) do not match!")
}

## initialize SGI structure; minsize is set to 5% of sample size
sg = sgi_init(hc, minsize = round(max(hc$order)/20), outcomes = selected_phenotypes) 
summary(sg) # shows how many samples in which cluster

## run SGI
as = sgi_run(sg)
as #prints significant associations, default adjusted p-value threshold is 0.05
#as$results
#summary(as, padj_th = 0.05)

# print(as, padj_th = 0.001) # different p-value
# print(as, padj_th = 0.1, by_clusters = F) #print with respect to outcomes 


# -----------------------------
# Find Samples in Cluster
# -----------------------------
cluster_assignments <- cutree(hc, h = 47.81320 )  # Adjust 'k' or 'h' based on significant clusters/heights of significant cluster pairs
table(cluster_assignments) #cluster sizes

plot(hc,labels = F)
rect.hclust(hc, h=47.81320 , border = 2:6)
abline(h=47.81320 , col = "red")

cluster_assignment <- data.frame(
  Sample = names(cluster_assignments),
  Cluster = cluster_assignments
)
#print(cluster_assignment)

samples_by_cluster <- split(cluster_assignment$Sample, cluster_assignment$Cluster)
#print(samples_by_cluster)

#outcome_with_sample_column <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_NOTwinsorized.xlsx", sheet = "Selected Phenotype Data"))
#final_data <- merge(cluster_assignment, outcome_with_sample_column, by = "Sample")
#cluster_single_level_249 <- subset(final_data, Cluster %in% c(1,76)s) # 
#write_xlsx(cluster_level_9, "cluster_level_9.xlsx")

# -----------------------------
# Plots
# -----------------------------
## generate tree plot, show results for adjusted p-values <0.05
gg_tree = plot(as, padj_th = 0.05)
gg_tree

#pdf("gg_tree_overview_ward_wind_clean_adc_eucl.pdf", width = 8, height = 6)

## plot overview, including clinical data and metabolomics data matrix
plot_overview( gg_tree = gg_tree, as = as, 
               outcomes = selected_phenotypes,
               outcomes_nmax = length(selected_phenotypes),
               xdata    = data_scaled,
               data_cluster_cols = 1,
               data.color = colorRampPalette(c("blue", "white", "red"))(100),
               data_title = "metabolites", 
               )

#dev.off()

ggs = plot_outcomes(sg, as, padj_th = 0.05) ## show that ggs is a list of plots 


### Heatmap only
library(pheatmap)

# neue Farbskala mit SD3 (weil sonst alles rot)
my_breaks <- seq(-3, 3, length.out = 101)

pheatmap(
  t(data_scaled),
  cluster_cols = hc,
  #breaks = my_breaks,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap (hclust results)"
)

### Interactive Heatmap (https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html)
library(heatmaply)

capped_data <- as.data.frame(pmax(pmin(t(data_scaled), 3), -3))
##davor wars (t(data_scaled) bei heatmaply ohne breaks -> aber dann rote Tabelle wie bei SGI

heatmaply(
  capped_data,
  Colv  = as.dendrogram(hc),
  colors = colorRampPalette(c("blue", "white", "red"))(100),  # color palette as in SGI
  main = "Interactive Heatmap +/- 3SD (hclust results)",
  file = "heatmaply_plot_wind_clean_cadj_3SD.html" ##is it 2 or 4?
)
#browseURL("heatmaply_plot_windsorized.html")

# -----------------------------
# make Plots automatically
# -----------------------------
results <- capture.output(as)
results <- results[!results %in% c("", "SGI associations...", "--- Significant results ---")]

# parsing results
sgi_results <- data.frame(line = results) %>%
  separate(line, into = c("Comparison", "Variables"), sep = " : ") %>%
  separate_rows(Variables, sep = ",\\s*")

combined_plots <- list()
unique_comp <- unique(sgi_results$Comparison)

for (comparison in unique_comp) {
  variables <- trimws(sgi_results$Variables[sgi_results$Comparison == comparison]) # alle Variable aiswählen die zu aktueller Comparison gehören - GETRIMMT WEIL LEERZEICHEN GR
  if (length(variables) == 1){
    combined_plots[[comparison]] <- ggs[[variables]][[comparison]]
  } else {
    combined_plots[[comparison]] <- Reduce("+", lapply(variables, function(variables) ggs[[variables]][[comparison]]))
  }
}

# combines plots (und speichern pdf)
#pdf("sgi_plots_combined_selected_phenotypes.pdf", width = 8, height = 6)
combined_plots
#dev.off()

# -----------------------------
# Plots manually
# -----------------------------
#pdf("boxplots_ward_adjusted.pdf", width = 8, height = 6)

# look at subgroups
#ggs$Country$`2vs3`
#ggs$Country$`4vs5`+ ggs$blood_ldl$`4vs5`+ ggs$blood_chol$`4vs5`+ ggs$blood_trig$`4vs5`+ ggs$blood_apo_ai$`4vs5`
#ggs$Country$`8vs9`+ ggs$blood_ldl$`8vs9`+ ggs$blood_trig$`8vs9`
#ggs$Country$`14vs15`+ ggs$waist_who$`14vs15`+ ggs$blood_leptin$`14vs15`
#ggs$Country$`16vs17`+ ggs$blood_trig$`16vs17`+ ggs$tg_hdl_ratio$`16vs17`
#ggs$Country$`18vs19`+ ggs$blood_ldl$`18vs19`+ ggs$blood_chol$`18vs19`+ ggs$blood_trig$`18vs19`+ ggs$blood_apo_b$`18vs19`+ ggs$blood_apo_ai$`18vs19`+ ggs$homa_index$`18vs19` + ggs$tg_hdl_ratio$`18vs19`
#ggs$Country$`22vs23`+ ggs$blood_chol$`22vs23`+ ggs$blood_igf1$`22vs23`+ ggs$blood_glucose$`22vs23`
#ggs$Country$`28vs29`+ ggs$ufa_perc$`28vs29`+ ggs$bmicatiotf$`28vs29`+ ggs$blood_hdl$`28vs29`+ ggs$blood_igf1$`28vs29` + ggs$tg_hdl_ratio$`28vs29`

#dev.off()
# ==================================
# kmeans (consensus clustering)
# ==================================

# -----------------------------
# find optimal k #### https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
# -----------------------------
# Elbow method 
fviz_nbclust(data_scaled, kmeans, k.max=15, method = "wss", print.summary = T) +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(data_scaled, kmeans, k.max = 15, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# verbose = FALSE to hide computing progression.
set.seed(123)
custom_kmeans <- function(data_scaled, centers, nstart) {
  kmeans(data_scaled, centers = centers, nstart = 25, iter.max = 100)
}
fviz_nbclust(data_scaled, custom_kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method")

## Nb clust
library(NbClust)
NbClust(data_scaled, distance = "euclidean", method = "kmeans")


#### Consensus Clustering ####

#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

results_consensus <- ConsensusClusterPlus(as.matrix(t(data_scaled)), 
                                          maxK = 15, 
                                          reps = 1000, 
                                          pItem = 0.8,
                                          pFeature = 1,
                                          seed = 4242,
                                          clusterAlg = "km", 
                                          distance = "euclidean", 
                                          plot = "png",
)

icl = calcICL(results_consensus,title="ICL",plot="png")
icl[["clusterConsensus"]]

#k3_clusters <- data.frame(results_consensus[[3]]$consensusClass)     
#colnames(k3_clusters)[1] <- "Cluster"

# -----------------------------
# actual kmeans clustering
# -----------------------------
library(clusterSim)

smallest_wcss <- Inf
index <- 0

for (i in 2:16){
  set.seed(4242)
  kmeans_results <- kmeans(data_scaled, centers = i, nstart = 25)
  #kmeans_results$clustes
  
  kmeans_cluster <- as.data.frame(kmeans_results$cluster)
  
  # Davies-Bouldin-Index: Kompaktheit und separation, je kleiner der wert desto besser
  davies_bouldin_index <- index.DB(data_scaled, kmeans_results$cluster)
  
  # Gib den Davies-Bouldin-Index aus
  print(paste("Davies-Bouldin-Index für k =", i, ":", davies_bouldin_index[1]))
}

set.seed(4242)
kmeans10 <- kmeans(data_scaled, centers = 10, nstart = 25)
kmeans_cluster10 <- as.data.frame(kmeans10$cluster)
table(kmeans10$cluster) ## where cluster 2 is the most stable (bootstrapping)
which(kmeans10$cluster == 7)


saveRDS(kmeans10, "kmeans_10_results.RDS")

#### Bootstrapping
library(fpc)
k <- 7:14
res_loop_boot <- data.frame(
  k = k, mean_jac = NA)

for (i in k){
  print(i)
  res_boot <- clusterboot(data_scaled,
                          runs = 25, #probiert 25 verschiedene initial starts and keeps the best
                          krange = i,
                          clustermethod = kmeansCBI,#kmeans
                          B = 1000,            
                          bootmethod = "boot",
                          multipleboot = TRUE, # punkte können mehrfach gezogen werden
                          bscompare = TRUE,
                          count = FALSE, # zeigt nicht counting am screen
                          seed = 4242)
  print(res_boot)
} 

# res_boot$partition: Cluster-Zuweisung 
# res_boot$result$clustercenters: Die Cluster-Zentren
# res_boot$bootresult: Infos zur Stabilität

#Visualize Jaccard indices (heatmap)
#heatmap(res_boot$bootresult$jaccard, Rowv = NA, Colv = NA, scale = "none", col = heat.colors(256), main = "Jaccard Index Heatmap")

idx_cluster2_TRUE <- which(res_boot$partition == 6)
idx_cluster2_TRUE    # Datenzeilen aus Cluster 2

##### silhouette score #####
sil <- silhouette(kmeans10$cluster, dist(data_scaled))
summary(sil)
fviz_silhouette(sil)
avg_sil <- mean(sil[, 3])
avg_sil


##### connectivity ###### 
# for local connectivity (kompaktheit) https://www.rdocumentation.org/packages/clValid/versions/0.7/topics/clValid
library(clValid)
conn <- clValid(data_scaled,
                nClust = 3:15,
                clMethods = c("kmeans", "hierarchical"),
                validation = "internal",
                metric = "euclidean", 
                method = "average",
)

summary(conn)


#saveRDS(idx_cluster2, "stabilstes_Cluster_Bootstrap.RDS")

# ==================================
# GMM https://cran.r-project.org/web/packages/ClusterR/vignettes/the_clusterR_package.html
# ==================================
#library(ClusterR) ## also try mclust()?
#gmm = GMM(data_scaled, 10, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
 #         em_iter = 10, verbose = F)   

# predict centroids, covariance matrix and weights
#pr =  predict_GMM(data = data_scaled,
#                  CENTROIDS = gmm$centroids,
 #                 COVARIANCE = gmm$covariance_matrices,
  #                WEIGHTS = gmm$weights)
#print(pr)
#cluster_assignments <- pr$cluster_labels
#cluster_assignments


library(mclust)
# chooses optimal cluster (G)
fit <- Mclust(data_scaled, G=1:14)
summary(fit)

fit$classification

plot(fit, what = c("classification"))
plot(fit, what = "BIC")


# speichern von cluster
#gmm_cluster <- data.frame(
#  cluster = cluster_assignments
#)
#rownames(gmm_cluster) <- rownames(data_scaled)

# ==================================
# Subspace Clustering
# ==================================
options(java.parameters = "-Xmx4g")
library(subspace)
library(rJava)

#xi = in wieviele Intervalle wird jeder metabolit unterteilt
#tau = "anzahl an datenpunkten damit Zelle als dicht gilt" 40% bei 0.4


## xi = 15, tau = 0.1: 762 Cluster max 3D !!!!
## xi = 15, tau = 0.2: 133 cluster max 2D
## xi = 15, tau = 0.3: 9 Cluster max 1-D
## xi = 15, tau = 0.4: 2 Cluster max 1D
## xi = 10, tau = 0.1: 14378 Cluster max 4D !!!!
## xi = 10, tau = 0.2: 470 Cluster max 2D
## xi = 10, tau = 0.3: 106 Cluster max 2D
## xi = 10, tau = 0.4: 10 Cluster max 1D
## xi = 9, tau = 0.2: 878 Cluster max 3D !!!!
## xi = 9, tau = 0.25: 392 Cluster max 3D !!!!
## xi = 9, tau = 0.3: 187 Cluster max 2D
## xi = 9, tau = 0.35: 72 Cluster max 2D
## xi = 9, tau = 0.4: 26 Cluster max 2D
## xi = 9, tau = 0.45: 11 Cluster max 2D
## xi = 9, tau = 0.5: 5 Cluster max 1D
## xi = 8, tau = 0.2: 1863 Cluster max 4D !!!!
## xi = 8, tau = 0.25: 675 Cluster max 3D !!!!
## xi = 8, tau = 0.3: 310 Cluster max 3D !!!!
## xi = 8, tau = 0.35: 146 Cluster max 2D
## xi = 8, tau = 0.4: 60 Cluster max 2D
## xi = 8, tau = 0.45: 18 Cluster max 2D
## xi = 8, tau = 0.5: 8 Cluster max 2D
## xi = 7, tau = 0.2: 4972 Cluster max 4D !!!! 1 Minute
## xi = 7, tau = 0.25: 1164 Cluster max 3D !!!!
## xi = 7, tau = 0.3: 407 Cluster max 3D !!!!
## xi = 7, tau = 0.35: 206 Cluster max 2D
## xi = 7, tau = 0.4: 95 Cluster max 2D
## xi = 7, tau = 0.45: 41 Cluster max 2D
## xi = 7, tau = 0.5: 15 Cluster max 1D
## xi = 6, tau = 0.2: 23409 Cluster max 4D !!!! 12 Minuten
## xi = 6, tau = 0.25: 5725 Cluster max 4D !!!!
## xi = 6, tau = 0.3: 1378 Cluster max 3D !!!!
## xi = 6, tau = 0.35: 443 Cluster max 2D
## xi = 6, tau = 0.4: 206 Cluster max 2D
## xi = 6, tau = 0.45: 103 Cluster max 2D
## xi = 6, tau = 0.5: 50 Cluster max 1D
## xi = 3, tau = 0.5; 28153 objects


system.time({
  clique_result <- CLIQUE(data_scaled, xi = 15, tau = 0.2)
})
clique_result
str(clique_result)
clique_result <- clique_result[order(unlist(lapply(lapply(clique_result, "[[", "subspace"),"sum")),decreasing = TRUE)] 

clusters_clique <- clique_result$clusters.subspace
clusters_clique

## loop?
xi_values <- c(20,30,40
               #,5,4,3
)     
tau_values <- c(0.1, 0.05, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5
)  


results_list <- list()
res_counter <- 0

for (xi in xi_values) {
  for (tau in tau_values) {
    cat("\n---------------------------------------\n")
    cat("Starte CLIQUE mit xi =", xi, ", tau =", tau, "...\n")
    
    start_time <- Sys.time()
    
    
    # CLIQUE mit Zeitlimit
    tryCatch({
      clique_res <- CLIQUE(data_scaled, xi = xi, tau = tau)
      
      end_time <- Sys.time()
      time_diff <- round(difftime(end_time, start_time, units = "mins"), 2)
      
      cat("\n---------------------------------------\n")
      cat("CLIQUE mit xi =", xi, ", tau =", tau, "beendet nach", time_diff, "Minuten.\n")
      
      # Überprüfe die Cluster und speichere die Ergebnisse (wie zuvor)
      dims_per_cluster <- sapply(clique_res, function(cl) sum(cl$subspace))
      max_dims <- max(dims_per_cluster)
      
      cat("Gefundene Cluster:", length(clique_res), 
          "--> maximale Anzahl Dimensionen in einem Cluster:", max_dims, "\n")
      
      if (any(dims_per_cluster >= 3)) {
        cat("Mindestens ein Cluster mit >= 3 Dimensionen gefunden! Speichere Ergebnis...\n")
        
        res_counter <- res_counter + 1
        results_list[[res_counter]] <- list(
          xi = xi,
          tau = tau,
          clique_result = clique_res,
          dims_per_cluster = dims_per_cluster
        )
      } else {
        cat("Kein Cluster mit >= 3 Dimensionen. NICHT abgespeichert.\n")
      }
      
    }, error = function(e) {
      end_time <- Sys.time()
      time_diff <- round(difftime(end_time, start_time, units = "mins"), 2)
      cat("\n---------------------------------------\n")
      cat("Fehler bei CLIQUE mit xi =", xi, ", tau =", tau, "nach", time_diff, "Minuten:", e$message, "\n")
      
    }, interrupt = function(i){
      end_time <- Sys.time()
      time_diff <- round(difftime(end_time, start_time, units = "mins"), 2)
      cat("\n---------------------------------------\n")
      cat("Subgrouping bei xi =", xi, "und tau =", tau, "dauert länger als 10 Minuten. (nach ",time_diff," Minuten) --> Abbruch\n")
    }
    )
  }
}


cat("\n\n--- Zusammenfassung ---\n")
cat("Insgesamt", res_counter, "Parameterkombinationen erzeugten mindestens einen 3D-Cluster.\n")

length(results_list)


######results

visualization <- list()

# loop through results
for (i in 1:length(results_list)) {
  result <- results_list[[i]]
  xi <- result$xi
  tau <- result$tau
  clique_result <- result$clique_result
  dims_per_cluster <- result$dims_per_cluster
  
  cat("\n---------------------------------------\n")
  cat("Ergebnisse für xi =", xi, ", tau =", tau, ":\n")
  
  # loop through clusters
  for (j in 1:length(clique_result)) {
    cluster <- clique_result[[j]]
    dimensionen <- which(cluster$subspace) # Dimensionen (Metaboliten), die TRUE sind
    anzahl_dimensionen <- length(dimensionen)
    
    # only consider clusters with more than 4 dimensions
    if (anzahl_dimensionen >= 4) {
      samples <- cluster$objects
      sample_names <- rownames(data_scaled)[samples]
      
      cat("  Cluster", j, ":\n")
      cat("    Dimensionen (Metaboliten):", paste(colnames(data_scaled)[dimensionen], collapse = ", "), "\n")
      cat("    Anzahl der Samples:", length(samples), "\n")
      #cat("    Sample-Indices:", paste(samples, collapse = ", "), "\n")
      cat("    Sample-Namen:", paste(rownames(data_scaled)[samples], collapse = ", "), "\n")
      
      # save data for visualisation
      visualization[[length(visualization) + 1]] <- list(
        xi = xi,
        tau = tau,
        cluster = j,
        dimensionen = colnames(data_scaled)[dimensionen],
        anzahl_samples = length(samples),
        sample_names = sample_names,
        metaboliten_werte = colMeans(data_scaled[samples, ])
      )
    }
  }
}


#saveRDS(results_list, file = "results_list_6_bis_9.rds")
rl <- readRDS("results_list_6_bis_9.rds")


# ==================================
# SOM (self organizing map) https://www.r-bloggers.com/2021/04/self-organizing-maps-in-r-supervised-vs-unsupervised/
# ==================================
library(kohonen)

# Grid definieren mit seed
set.seed(42)
grid <- somgrid(xdim = 4, ydim = 4, topo = "hex") # xdim und ydim größe des gitters, topo= sechseckig oder rechteckig ("rect")

# SOM trainieren
data_matrix <- as.matrix(data_scaled)
map <- som(
  data_matrix,
  grid = grid,
  #rlen = 100, # = epochs
  #alpha = c(0.05, 0.01),# = lernrate (start und endwerte)
  #keep.data = TRUE #damit plot funktionen genutzt werden könne (?)
)

# code map (?)
plot(map)

# Hitmap
plot(map, type = "counts")

# Classification
zuordnungen <- map$unit.classif
print(zuordnungen)

# show per neuron
daten_pro_neuron <- split(data_scaled, zuordnungen)

# Zugriff auf Datenpunkte für ein bestimmtes Neuron (z. B. Neuron 5)
daten_neuron_5 <- daten_pro_neuron[["5"]]

# distance map
plot(map, type = "dist.neighbours") #Visualisiert die Distanzen zwischen benachbarten Neuronen, Hohe Werte (dunkle Farben) in der U-Matrix deuten auf Clustergrenzen hin.

#evaluierung
# Quantisierungsfehler
quant_fehler <- sum(map$distances) / nrow(data_scaled)

# Topologischer Fehler
topo_fehler <- topographic.error(map)
