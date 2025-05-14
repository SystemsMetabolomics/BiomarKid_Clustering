source("00_init.R")
source("useful_functions.R")
source("utils_plotting.R")
source("utils_data_preparation.R")

# ==================================
# Outlier Detection with PCA and mahalanobis distance
# ==================================
## calculate outliers with mahalanobis distance (across all compartments of pca) 
outlier_mahalanobis <- function(pca_scores, p_cutoff) {
  mahal_dist <- mahalanobis(pca_scores, colMeans(pca_scores), cov(pca_scores))
  cutoff <- qchisq(p_cutoff, ncol(pca_scores))
  outliers <- which(mahal_dist > cutoff)
  return (outliers)
} 

outlier_detection <- function (data,
                               phenotypes,
                               chi_cutoff = 0.9999, 
                               save_outlier_plot = TRUE,
                               create_exp_plots = FALSE) {
  
  data_log_z <- data_preprocessing(data)
  data_adj <- adjust_for_country(data_log_z, phenotypes)
  
  pca <- prcomp(data_adj, scale. = F) # scale = false because data is already scaled
  scores <- pca$x  # PCA-transformed data
  
  outliers <- outlier_mahalanobis(scores, chi_cutoff)
  
  if (length(outliers) > 0){
  if (save_outlier_plot) {
    plot_pca_outliers(scores,
                      outliers,
                      paste0("output/plots/outlier_detection_", chi_cutoff, ".pdf")) 
  }
  
  ####### Data Export #######
  #### remove outliers from raw data files and create new ones
  data_no_outliers <- data[-outliers, ]
  pheno_no_outliers <- phenotypes[!rownames(phenotypes) %in% outliers, ]
  
  saveRDS(data_no_outliers, "data/raw/metabolites_no_outliers.rds")
  saveRDS(pheno_no_outliers,"data/raw/phenotypes_no_outliers.rds")
  
  ####### Exploratory Plots #######
  if (create_exp_plots) {
    #### for PCA
    ## scree plot
    p1 <- fviz_eig(pca, addlabels = TRUE) # only 30% in first 2 Components
    print (p1)
    ggsave(filename = "output/plots/scree_plot_PCA.png",
           plot = p1)
    
    # BiPlot (https://www.datacamp.com/de/tutorial/pca-analysis-r)
    print(fviz_pca_var(pca, col.var = "black")) # metabolites/Pfeile
    print(fviz_pca_ind(pca, col.ind = "black", labels = T)) # dots only (samples)
  }
  return (list(data_no_outliers, pheno_no_outliers))
  } else {
    message("No outliers found. Turn off Outlier Detection.")
  }
}
