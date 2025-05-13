source("00_init.R")
source("utils_phenotype_association.R")
source("useful_functions.R")
source("utils_T132_phenotype_selection.R")

cluster_hc <- readRDS("output/cluster_hc_list.RDS")
cluster_km <- readRDS("output/cluster_km_list_nstart150.RDS")
cluster_gmm <- readRDS("output/cluster_gmm_list.RDS")
cluster_som <- readRDS("output/cluster_som_list.RDS")
data_scaled <- readRDS("data/processed/metabolites_clean_no_outliers.rds")


# ==================================
# CHOOSE SETTINGS
# ==================================
##choose timepoint of phenotypes for phenotype associationoca
timepoint = "T132" #(T66 or T132)

## select k_values phenotype associations should be done on (later different sheets in excel file)
k_values <- c(7, 13)

## select final K (or k plots should be done with -> needs to be part of k_values
k <- 13

# should SGI toolbox be run and create plots as pdf?
run_SGI <- TRUE
# ==================================
#
# ==================================

if (timepoint == "T66"){
  # load phenotypes same timepoint
  outcomes_T66 <- readRDS("data/processed/selected_phenotypes_no_outliers.rds")
  outcomes_T66$Sample <- rmv_tp(outcomes_T66$Sample)
  outcomes_T66 <- samples_as_rownames(outcomes_T66)
  outcomes_log <- outcomes_T66
} else if (timepoint == "T132"){
  # load phenotypes of both timepoints 
  later_phenotype <- readRDS("data/raw/phenotypes_T132.rds")
  T66 <- readRDS("data/processed/selected_phenotypes_no_outliers.rds")
  create_selected_phenotypes_for_later_timepoint(T66, later_phenotype)
  outcomes_T132 <- readRDS("data/processed/selected_phenotypes_T132_as_in_T66.rds")
  outcomes_T132 <- samples_as_rownames(outcomes_T132)
  outcomes_log <- outcomes_T132
}

# log2 tranformation of outcome values
outcomes_log <- outcomes_log %>%
  mutate(across(where(is.numeric), ~log2(.)))

outcomes_cluster_km <- merge_cluster_assignments(outcomes_log, cluster_km, k_values)
outcomes_cluster_gmm <- merge_cluster_assignments(outcomes_log, cluster_gmm, k_values)
outcomes_cluster_hc <- merge_cluster_assignments(outcomes_log, cluster_hc, k_values)
outcomes_cluster_som <- merge_cluster_assignments(outcomes_log,cluster_som, k_values)


# Welch's T-Test from utils_phenotype_association file
pairwise_tests_km <- cluster_test(outcomes_log, outcomes_cluster_km, k_values, compare_type = "pairwise" )
one_vs_all_tests_km <- cluster_test(outcomes_log, outcomes_cluster_km, k_values, compare_type = "one_vs_all" )
pairwise_tests_hc <- cluster_test(outcomes_log, outcomes_cluster_hc, k_values, compare_type = "pairwise" )
one_vs_all_tests_hc <- cluster_test(outcomes_log, outcomes_cluster_hc, k_values, compare_type = "one_vs_all" )
pairwise_tests_gmm <- cluster_test(outcomes_log, outcomes_cluster_gmm, k_values, compare_type = "pairwise" )
one_vs_all_tests_gmm <- cluster_test(outcomes_log, outcomes_cluster_gmm, k_values, compare_type = "one_vs_all" )
pairwise_tests_som <- cluster_test(outcomes_log, outcomes_cluster_som, k_values, compare_type = "pairwise" )
one_vs_all_tests_som <- cluster_test(outcomes_log, outcomes_cluster_som, k_values, compare_type = "one_vs_all" )


list_tests <- list(pairwise_km = pairwise_tests_km,
                      pairwise_hc = pairwise_tests_hc,
                      paiwise_gmm = pairwise_tests_gmm,
                      paiwise_som = pairwise_tests_som,
                      one_vs_all_km = one_vs_all_tests_km,
                      one_vs_all_hc = one_vs_all_tests_hc,
                      one_vs_all_gmm = one_vs_all_tests_gmm,
                      one_vs_all_som = one_vs_all_tests_som)
                      
                      
save_test_results(list_tests)


# ==================================
# PLOTS DONT WORK FOR DUMMY DATA BECAUSE THERE ARE TOO MANY METHODS WITH NO SIGNIFICANCE
# ==================================
### PLOTS FOR FINAL CHOSEN K (or K one wants to look at - must be part of k_values)
if (k %in% k_values){
  # one vs all plots for only looked at cluster vs rest
  boxplot_one_vs_all_hc <- boxplots_sig_one_vs_all(outcomes_cluster_hc, outcomes_log, one_vs_all_tests_hc, "hc", k, timepoint= timepoint, TRUE)
  boxplot_one_vs_all_km <- boxplots_sig_one_vs_all(outcomes_cluster_km, outcomes_log, one_vs_all_tests_km, "km", k,timepoint= timepoint, TRUE)
  boxplot_one_vs_all_som <- boxplots_sig_one_vs_all(outcomes_cluster_som, outcomes_log, one_vs_all_tests_som, "som", k,timepoint= timepoint, TRUE)

  # boxplots for all cluster - one vs all is highlightes, for pairwise sig lines turn sig_lines to TRUE
  boxplot_km_13 <- boxplots_sig(outcomes_cluster_km, outcomes_log, "km", k,  timepoint= timepoint,save_plots = TRUE, sig_lines = FALSE)
  boxplot_som_13 <- boxplots_sig(outcomes_cluster_som, outcomes_log, "som", k,  timepoint= timepoint,save_plots = TRUE, sig_lines = FALSE)
  boxplot_hc_13 <- boxplots_sig(outcomes_cluster_hc, outcomes_log, "hc", k,  timepoint= timepoint,save_plots = TRUE, sig_lines = FALSE)

  # distribution plots for all cluster - one vs all is highlightes, for pairwise sig lines turn sig_lines to TRUE
  distribution_km_13 <- distributionplots_sig(outcomes_cluster_km, outcomes_log, "km", k,  timepoint= timepoint,TRUE, sig_lines = FALSE)
  distribution_hc_13 <- distributionplots_sig(outcomes_cluster_hc, outcomes_log, "hc", k,  timepoint= timepoint,TRUE,  sig_lines = FALSE)
  distribution_som_13 <- distributionplots_sig(outcomes_cluster_som, outcomes_log, "som", k,  timepoint= timepoint,TRUE, sig_lines = FALSE)
} else {
  message("Plots cannot be created because t-tests were not executed (k is not part of k_range).")
}


# ==================================
# SGI
# ==================================
if (run_SGI){
  if (timepoint == "T66"){
   selected_phenotypes <- outcomes_T66
  } else if (timepoint == "T132"){
    selected_phenotypes <- outcomes_T132
  }
  
  data_scaled_no_timepoint <- rownames_as_column(data_scaled)
  data_scaled_no_timepoint$Sample <- rmv_tp(data_scaled_no_timepoint$Sample)
  data_scaled_no_timepoint <- samples_as_rownames(data_scaled_no_timepoint)
  hc_ward2 = hclust(dist(data_scaled_no_timepoint), method = "ward.D2")
  hc_sgi <- hc_ward2
  
  ## refactoring values that are char
  selected_phenotypes <- selected_phenotypes %>%
    mutate(across(where(is.character), as.factor))
  
  if (!isTRUE(all.equal(hc_sgi$labels, rownames(selected_phenotypes)))) {
    stop("Identifiers in hc_sgi$labels and rownames(selected_phenotypes) do not match!")
  }
  
  ## initialize SGI structure; minsize is set to 5% of sample size
  sg = sgi_init(hc_sgi, minsize = round(max(hc_sgi$order)/20), outcomes = selected_phenotypes) 
  summary(sg) # shows how many samples in which cluster
  
  ## run SGI
  as = sgi_run(sg)
  as #prints significant associations, default adjusted p-value threshold is 0.05
  
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
  
  # combines plots (and safe pdf)
  pdf("output/plots/sgi_plots_combined_selected_phenotypes_wardD2.pdf", width = 8, height = 6)
  combined_plots
  dev.off()
}




