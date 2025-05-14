source("00_init.R")
source("utils_plotting.R")
source("utils_metabolite_relevance.R")

# ==================================
# CHOOSE SETTINGS
# ==================================

k = 13 # final clusters for heatmaps

# ==================================
# 
# ==================================

data_scaled <- readRDS("data/processed/metabolites_clean_no_outliers.rds")
cluster_hc <- readRDS("output/cluster_hc_list.RDS")
cluster_km <- readRDS("output/cluster_km_list_nstart150.RDS")
cluster_gmm <- readRDS("output/cluster_gmm_list.RDS")
cluster_som <- readRDS("output/cluster_som_list.RDS")

# ==================================
# Heatmaply for all samples and metabolites (HTML)
# ==================================
## color palette only for 12 colors (12 is okay - grey is added)
heatmap_every_sample_by_rows(data_scaled, cluster_km, method = "kmeans", k = k)
heatmap_every_sample_by_rows(data_scaled, cluster_gmm, method = "gmm",  k = k)
heatmap_every_sample_by_rows(data_scaled, cluster_hc, method = "hc", k = k)
heatmap_every_sample_by_rows(data_scaled, cluster_som, method = "som", k = k)

# ==================================
# Circlize Heatmap for cluster means
# ==================================
cluster_means_hc <- create_cluster_means(data_scaled, cluster_hc, k)
cluster_means_som <- create_cluster_means(data_scaled, cluster_som, k)
cluster_means_km <- create_cluster_means(data_scaled, cluster_km, k) # same as cluster centers from kmeans function

## not used in thesis because it's hard to present in PDF - circlize instead
# heatmaply(t(cluster_means_km),
#                colors = colorRampPalette(c("red", "white", "blue")),
#                Colv = FALSE,
#                show_dendrogram = c(FALSE,FALSE),
#                main = "Heatmap of Scaled Cluster Means (K-Means Clustering)",
#                fontsize_row = 10,
#                fontsize_col = 12,
#                height = 4000,
#                 #width = 1000,
#                scale = "row",
#                margins = c(80,160),
#                file = "Cluster_means_heatmap_km.html"
# )

# cluster metabolites with hc
row_hclust <- hclust(dist(t(data_scaled)), method = "ward.D2")
row_order <- rownames(t(data_scaled))[row_hclust$order]

pdf("output/plots/circos_heatmap_km.pdf")
plot_cluster_circos(cluster_means_km, row_order)
dev.off()

pdf("output/plots/circos_heatmap_hc.pdf")
plot_cluster_circos(cluster_means_hc, row_order)
dev.off()

pdf("output/plots/circos_heatmap_som.pdf")
plot_cluster_circos(cluster_means_som, row_order)
dev.off()


# ==================================
# Ranking
# ==================================
rank_metabolites(cluster_means_km, top=10, bottom = 10, paste0("output/ranked_metabolites_per_cluster", k, "_km.csv"))
rank_metabolites(cluster_means_hc, top=10, bottom = 10, paste0("output/ranked_metabolites_per_cluster", k, "_hc.csv"))
rank_metabolites(cluster_means_som, top=10, bottom = 10, paste0("output/ranked_metabolites_per_cluster", k, "_som.csv"))


