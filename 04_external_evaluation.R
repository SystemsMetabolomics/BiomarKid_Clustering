source("00_init.R")
source("utils_plotting.R")

# ==================================
# CHOOSE SETTINGS
# ==================================
k_range <- 2:16
k = 13 # chosen k for final clustering
# ==================================
# 
# ==================================

evaluation_hc <- readRDS("output/evaluation_hc_results.RDS")
evaluation_km <- readRDS("output/evaluation_km_results_nstart150.RDS")
evaluation_gmm <- readRDS("output/evaluation_gmm_results.RDS")
evaluation_som <- readRDS("output/evaluation_som_results.RDS")

cluster_hc <- readRDS("output/cluster_hc_list.RDS")
cluster_km <- readRDS("output/cluster_km_list_nstart150.RDS")
cluster_gmm <- readRDS("output/cluster_gmm_list.RDS")
cluster_som <- readRDS("output/cluster_som_list.RDS")

methods <- list("hc" = cluster_hc, "km" = cluster_km, "gmm" = cluster_gmm, "som" = cluster_som)

# ==================================
# Evaluation Plots 4 methods
# ==================================
## creates evaluation plots of 4 methods to compare and chose k
evaluations <- list(evaluation_hc, evaluation_km, evaluation_gmm, evaluation_som)
names(evaluations) <- names(methods)

pdf("output/plots/evaluation_plots.pdf", width = 8, height = 6)
evaluation_plots(evaluations, names(methods), "clustering methods")
dev.off()

# ==================================
# ARI Plots 4 methods
# ==================================
ari_plot <- ari_plots_methods(methods)
ggsave(paste0("output/plots/ari_plot.png"), width =8, height = 6, plot = ari_plot)

# ==================================
# Alluvial Plots
# ==================================
## add methods of clustering as prefix to clusternumber for alluvial plot to work
cluster_hc[] <- lapply(cluster_hc, function(x) paste0("hc", x))
cluster_km[] <- lapply(cluster_km, function(x) paste0("km", x))
cluster_gmm[] <- lapply(cluster_gmm, function(x) paste0("gmm", x))
cluster_som[] <- lapply(cluster_som, function(x) paste0("som", x))

p <- alluvial_plot(cluster_hc, cluster_km, cluster_gmm, cluster_som, k)
ggsave(paste0("output/plots/Alluvial_plot_methods_k_som", k, ".png"), width =8, height = 6, plot = p)

