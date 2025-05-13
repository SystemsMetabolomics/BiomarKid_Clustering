# weise cluster den daten zu
merge_cluster_assignments <- function(outcomes, cluster_assignments, cluster_cols){
  cluster_cols <- paste0("k_", cluster_cols)
  
  # remove timepoint from Sample ID
  rownames(cluster_assignments) <- rmv_tp(rownames(cluster_assignments))
  
  # create subset with only chosen k (cluster_cols)
  subset <- cluster_assignments[, cluster_cols, drop = FALSE]
 
  
  # change cluster assignment from numeric to factor
  subset <- mutate_if(subset, is.numeric, as.factor)
  
  # create Sample ID column to later merge by that column
  subset <- rownames_as_column(subset)
  outcomes <- rownames_as_column(outcomes)

  join_left <- left_join(outcomes, subset, by = "Sample")
  join_left <- samples_as_rownames(join_left)
  
  return(join_left)
}

# ==================================
# Pairwise t-test with bonferroni/Benjamin Hoch Correction
# ==================================
# numerical 
cluster_ttest <- function(numerical_data, data_cluster, k, compare_type = c("pairwise", "one_vs_all")){
  compare_type = match.arg(compare_type)
  results_list <- list()
  for (col_name in colnames(numerical_data)){
    cluster <- paste0("k_", k)
    omit_NA <- na.omit(data.frame(phenotype = numerical_data[[col_name]], cluster = data_cluster[[cluster]]))
    
    results_df <- data.frame(
      cluster_A = integer(),
      cluster_B = integer(),
      p_value = numeric(),
      adjusted_p_value = numeric(),
      significance = character(),
      phenotype = character()
    )
    
    if (compare_type == "pairwise"){
      group_counts <- table(omit_NA$cluster)
      singleton_clusters <- names(group_counts[group_counts < 2])
      if (length(singleton_clusters) > 0){
        omit_NA <- omit_NA[!omit_NA$cluster %in% singleton_clusters, ]
        message("Some Clusters were excluded from pairwise t-test because they consisted of only one Sample.")
      }
      
      if(k-length(singleton_clusters) < 2){
        message("There is only one cluster that is not a singleton cluster. Pairwise Test is not possible.")
      } else{
      pairwise_test <- pairwise.t.test(omit_NA$phenotype, omit_NA$cluster, p.adjust.method = "BH", paired = FALSE, pool.sd = FALSE)
      p_table <- as.data.frame(pairwise_test$p.value)
      for (i in 1:nrow(p_table)) {
        for (j in 1:(i)) {
          # wegen Dreieck
          p_value <- p_table[i, j]
          results_df <- rbind(results_df, data.frame(
            cluster_A = i + 1,
            cluster_B = j,
            p_value = NA,
            adjusted_p_value = p_value,
            significance = ifelse(p_value < 0.05, "significant", "not significant"),
            phenotype = col_name
          ))
        }
      }
      }
    } else if (compare_type == "one_vs_all"){
      unique_clusters <- unique(omit_NA$cluster)
      for (i in unique_clusters){
        i_cluster <- omit_NA$phenotype[omit_NA$cluster == i]
        rest <- omit_NA$phenotype[omit_NA$cluster != i]
        if (length(i_cluster) > 1 & length(rest) > 1) {
          ttest <- t.test(i_cluster, rest, var.equal = FALSE) ## var.equal = FALSE for non pooled SD's (Welch's Test)
          results_df <- rbind(results_df, data.frame(
            cluster_A = i,
            cluster_B = "others",
            p_value = ttest$p.value,
            adjusted_p_value = NA,
            significance = "not_calculated",
            phenotype = col_name
          ))
        } else {
          message(paste("Cluster", i, "has less than 2 values, t-Test not possible."))
        }
      }
      results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
      results_df$significance <- ifelse(results_df$adjusted_p_value < 0.05, "significant", "not_significant")
    }
    results_list[[col_name]] <- results_df
  }
  significant_results <- list()
  for (i in 1:length(results_list)) {
    filtered_values <- results_list[[i]] %>%
      filter(significance == "significant")
    if (nrow(filtered_values) >= 1) {
      significant_results <- rbind(significant_results, filtered_values)
    }
  }
  return(significant_results)
}

#categorical
cluster_cat_test <- function(categorial_data, data_cluster, k, compare_type = c("pairwise", "one_vs_all")){
  compare_type <- match.arg(compare_type)
  results_list <- list()
  for (pheno in colnames(categorial_data)){
    results_df <- data.frame(
      cluster_A = integer(),
      cluster_B = integer(),
      p_value = numeric(),
      phenotype = character(),
      test = character(),
      adjusted_p_value = numeric(),
      significance = character()
    )
    
    cluster <- paste0("k_", k)
    omit_NA <- na.omit(data.frame(phenotype = categorial_data[[pheno]], cluster = data_cluster[[cluster]]))
    
    unique_clusters <- unique(omit_NA$cluster)
    cont_table <- table(omit_NA$cluster, omit_NA$phenotype)
    
    if (compare_type == "pairwise") {
      k_clusters <- length(unique_clusters)
      
      for (i in 1:(k_clusters - 1)) {
        for (j in (i + 1):k_clusters) {
          pairwise_cont_table <- cont_table[c(i, j), ]
          
          
          suppressWarnings({ ## because I don't use these results anyways, just want to get the expected values
            chi_test <- chisq.test(pairwise_cont_table)
          })
          expected <- chi_test$expected
          
          test <- ifelse(any(expected < 5), "Fisher's Exact Test", "Chi-Square Test")
          test_result <- if (any(expected < 5)) fisher.test(pairwise_cont_table) else chisq.test(pairwise_cont_table)
          
          results_df <- rbind(
            results_df,
            data.frame(
              cluster_A = i,
              cluster_B = j,
              p_value = test_result$p.value,
              phenotype = pheno,
              test = test,
              adjusted_p_value = 0,
              significance = "not calculated"
            )
          )
        }
      }
    } else if (compare_type == "one_vs_all"){
      for (i in unique_clusters) {
        tmp_df <- omit_NA
        tmp_df$cluster <- ifelse(tmp_df$cluster == i, as.character(i), "others")
        
        cont_table_final <- table(tmp_df$cluster, tmp_df$phenotype)
        
        suppressWarnings({ ## because I don't use these results anyways, just want to get the expected values
          chi_test <- chisq.test(cont_table_final)
        })
        expected <- chi_test$expected
        
        test <- ifelse(any(expected < 5), "Fisher's Exact Test", "Chi-Square Test")
        test_result <- if (any(expected < 5)) fisher.test(cont_table_final) else chisq.test(cont_table_final)
        
        results_df <- rbind(
          results_df,
          data.frame(
            cluster_A = as.character(i),
            cluster_B = "others",
            p_value = test_result$p.value,
            phenotype = pheno,
            test = test,
            adjusted_p_value = 0,
            significance = "not calculated"
          )
        )
      }
    }
    results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
    results_df$significance <- ifelse(results_df$adjusted_p_value < 0.05, "significant", "not_significant")
    results_list[[pheno]] <- results_df
  }
  
  pairwise_significant_cat_results <- data.frame(
    cluster_A = character(),
    cluster_B = character(),
    p_value = numeric(),
    phenotype = character(),
    test = character(),
    adjusted_p_value = numeric(),
    significance = character()
  )
  
  for (pheno in names(results_list)) {
    filtered_values <- results_list[[pheno]] %>%
      filter(significance == "significant")
    if (nrow(filtered_values) >= 1) {
      pairwise_significant_cat_results <- rbind(pairwise_significant_cat_results, filtered_values)
    }
  }
  return(pairwise_significant_cat_results)
}

# decide which test
cluster_test <- function(data, data_cluster, k_values, compare_type = c("pairwise", "one_vs_all")) {
  compare_type <- match.arg(compare_type)
  numerical_data <- select_if(data, is.numeric)
  categorical_data <- select_if(data, is.factor)
  all_results <- list()
  
  for (k in k_values) {
    num_results <- list()
    cat_results <- list()
    
    if (ncol(numerical_data) > 0) {
      num_results <- cluster_ttest(numerical_data, data_cluster, k, compare_type)
    }
    
    if (ncol(categorical_data) > 0) {
      cat_results <- cluster_cat_test(categorical_data, data_cluster, k, compare_type)
    }
  
    ### only merge cluster if both have results, otherwise just state one or none with message
    n_num <- if (inherits(num_results, "data.frame")) nrow(num_results) else 0
    n_cat <- if (inherits(cat_results, "data.frame"))  nrow(cat_results)  else 0

    if (n_num > 0 && n_cat > 0){
      all_results[[as.character(k)]] <- bind_rows(num_results, cat_results)
    } else if (n_num > 0){
      all_results[[as.character(k)]] <- num_results
    } else if (n_cat > 0){
      all_results[[as.character(k)]] <- cat_results
    } else {
      message(paste("No significant clusters have been found for k =", k))
    }
  }
  return(all_results)
}

#### save t-test results as excel sheets
save_test_results <- function(list_tests){
  for (method in names(list_tests)){
    save_results <- list_tests[[method]]
    wb <- createWorkbook()
    # each dataframe its own sheet
    for (i in seq_along(save_results)) {
      addWorksheet(wb, sheetName = names(save_results)[i])
      writeData(wb, sheet = i, x = save_results[[i]])
    }
    file_name <- paste0("output/", method, "_", timepoint, ".xlsx")
    saveWorkbook(wb, file_name, overwrite = TRUE)
  }
}


# ==================================
# PLOTS
# ==================================
# adds lines and annotation if pairwise clusters are significantly different mit geomsignif
add_significant_lines <- function(significant_pairs, plot, y_max = NULL){
  comparisons_list <- lapply(1:nrow(significant_pairs), function(i) {
    c(
      as.character(significant_pairs$cluster_A[i]),
      as.character(significant_pairs$cluster_B[i])
    )
  })
  
  annotations_list <- sapply(significant_pairs$adjusted_p_value, function(p) {
    if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p < 0.05) {
      return("*")
    } else {
      return("ns")
    }
  })
  
  if (is.null(y_max)){
    plot <- plot +
      geom_signif(
        comparisons = comparisons_list,
        annotations = annotations_list,
        step_increase = 0.1
      )+
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = "*: p < 0.05\n**: p < 0.01\n***: p < 0.001",
        hjust = 1.1,
        vjust = 1.1,
        size = 3,
        color = "black"
      )} else {
        plot <- plot +
          geom_signif(
            comparisons = comparisons_list,
            annotations = annotations_list,
            step_increase = 0.25,
            aes(y = y_max * 0.9),
            tip_length = 0.03
          ) +
          annotate(
            "text",
            x = Inf,
            y = y_max * 1.4,
            label = "*: p < 0.05\n**: p < 0.01\n***: p < 0.001",
            hjust = 1.1,
            vjust = 1.1,
            size = 3,
            color = "black"
          )
      } 
  print(plot)
}

boxplots_sig_one_vs_all <- function(outcome_cluster, phenotypes, one_vs_all, clustering_method, k, timepoint, save_plots = FALSE){
  phenotypes_num <- select_if(phenotypes, is.numeric)
  k_cluster_str <- paste0("k_", k)
  
  one_vs_all_data <- one_vs_all[[as.character(k)]]
  num_sig <- one_vs_all_data[is.na(one_vs_all_data$test),]
  
  if (save_plots){pdf(paste0("output/plots/boxplots_one_vs_all_", clustering_method, "_", k, "_", timepoint, ".pdf"))}
  
  for (i in seq_len(nrow(num_sig))){
    pheno <- num_sig$phenotype[i]
    cluster_id <- num_sig$cluster_A[i]
    
    # chose current phenotype and k
    plot_data <- outcome_cluster[, c(pheno, k_cluster_str)]
    colnames(plot_data) <- c("value", "cluster")
    
    # Label the cluster vs all other clusters
    plot_data$group <- ifelse(plot_data$cluster == cluster_id, paste0("Cluster_", cluster_id), "all others")
    plot_data <- plot_data[!is.na(plot_data$value), ]
    
    fill_colors <- c("darkgreen", "black")
    names(fill_colors) <- c(paste0("Cluster_", cluster_id), "all others")
    
    plot <- ggplot(plot_data, aes(x = group, y = value, fill = group, color = group)) +
      geom_boxplot(width = 0.3, alpha = 0.4) +
      geom_jitter(width = 0.1, alpha = 0.8) +
      scale_fill_manual(values = fill_colors) +
      scale_color_manual(values = fill_colors) +
      labs(
        title = paste0("Boxplot: ", pheno, " (Cluster ", cluster_id, " vs. Others)"),
        x = NULL,
        y = pheno
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    print(plot)
  }
  
  if (save_plots) {
    dev.off()
  }
}

# boxplots for numerical phenotypes
boxplots_sig <- function(outcome_cluster, phenotypes, clustering_method, k, timepoint, save_plots = FALSE, sig_lines = TRUE) {
  phenotypes_num  <- select_if(phenotypes, is.numeric)
  k_cluster_str <- paste0("k_", k)
  
  if (save_plots){pdf(paste0("output/plots/boxplots_", clustering_method, "_", k, "_", timepoint, ".pdf"))}
  for (pheno in colnames(phenotypes_num)) {
    outcome_cluster_na <- outcome_cluster[!is.na(outcome_cluster[, pheno]), ]
    
    one_vs_all_tests <- paste0("one_vs_all_tests_", clustering_method)
    pairwise_tests <- paste0("pairwise_tests_", clustering_method)
    
    one_vs_all_tests_data <- get(one_vs_all_tests)
    pairwise_tests_data <- get(pairwise_tests)
    
    significant_cluster <- one_vs_all_tests_data[[as.character(k)]] %>%
      filter(phenotype == pheno)
    
    significant_cluster_ids <- if (nrow(significant_cluster) > 0) {
      significant_cluster$cluster_A
    }
    
    significant_pairs <- pairwise_tests_data[[as.character(k)]] %>%
      filter(phenotype == pheno)
    
    plot <- ggplot(outcome_cluster_na, aes_string(x = paste0("factor(", k_cluster_str, ")"), y = pheno)) +
      geom_boxplot(aes_string(fill = paste0("factor(", k_cluster_str, " %in% significant_cluster_ids)")), width = 0.5, alpha = 0.4) +
      scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "transparent"), guide = "none") +
      geom_jitter(width = 0.1, alpha = 0.5) +
      labs(title = paste0("Boxplots for ", pheno, " (", clustering_method, ", k=", k, ")"), x = "cluster", y = pheno) +
      theme_minimal()
    
    if (nrow(significant_pairs) > 0 && sig_lines) {
      add_significant_lines(significant_pairs, plot)
    } 
      print(plot)
  }
  if (save_plots) {dev.off()}
}

# distribution plots for categorial phenotypes
distributionplots_sig <- function(outcome_cluster, phenotypes, clustering_method, k, timepoint, save_plots = TRUE, sig_lines = TRUE){
  phenotypes_cat <- select_if(phenotypes, is.factor)
  k_cluster_str <- paste0("k_", k)
  
  if (save_plots) {pdf(paste0("output/plots/Distribution_Plot_", clustering_method, "_", k, "_", timepoint, ".pdf"), width = 8, height = 8)}
  for (pheno_name in names(phenotypes_cat)) {
    # Calculate cluster sizes for y-axis scaling
    cluster_sizes <- rowSums(table(outcome_cluster[[k_cluster_str]], outcome_cluster[[pheno_name]]))
    max_y <- max(cluster_sizes)
    
    one_vs_all_tests <- paste0("one_vs_all_tests_", clustering_method)
    pairwise_tests <- paste0("pairwise_tests_", clustering_method)
    
    one_vs_all_tests_data <- get(one_vs_all_tests)
    pairwise_tests_data <- get(pairwise_tests)
    
    
    significant_cluster <- one_vs_all_tests_data[[as.character(k)]] %>%
      filter(phenotype == pheno_name)
    
    significant_cluster_ids <- if (nrow(significant_cluster) > 0) {
      significant_cluster$cluster_A
    }
    
    plot <- ggplot(outcome_cluster, aes(x = factor(!!sym(k_cluster_str)), fill = .data[[pheno_name]])) +
      geom_bar(aes(y = after_stat(count), colour = factor(!!sym(k_cluster_str)) %in% significant_cluster_ids), stat = "count", linewidth = 1.2) +
      scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "transparent"), guide = "none") +
      labs(
        title = paste("Distribution of", pheno_name, "in", k, "cluster"),
        x = "cluster",
        y = "number",
        fill = pheno_name
      ) +
      theme_bw()
    
    significant_pairs <- pairwise_tests_data[[as.character(k)]] %>%
      filter(phenotype == pheno_name)
    
    if (nrow(significant_pairs) > 0 && sig_lines) {
      add_significant_lines(significant_pairs, plot, max_y)
    } 
      print(plot)
  }
  if (save_plots) {dev.off()}
}