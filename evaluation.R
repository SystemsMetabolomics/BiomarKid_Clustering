source("data_preprocessing.R")

library(Rtsne)
library(cluster)

outcomes <- as.data.frame(read_excel("/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx", sheet = "Selected Phenotype Data"))
rownames(outcomes) <- outcomes[,1]
outcomes <- outcomes[,-1]

kmeans_results10 <- readRDS("/Users/denise/Desktop/BA/kmeans_10_results.RDS")
kmeans_results8 <- readRDS("/Users/denise/Desktop/BA/kmeans_8_results.RDS")

cluster_kmeans <- kmeans_results10$cluster

table(cluster_kmeans)


#### NUMERICAL
outcomes_num <- select_if(outcomes, is.numeric)
comb_outcomes_cluster10 <- outcomes_num
comb_outcomes_cluster10$cluster <- cluster_kmeans

####### Prüfen ob Phenotypes normalverteilt sind ######
## visuell mit histogramm
for (col_name in names(outcomes_num)) {
  if (is.numeric(outcomes[[col_name]])) {
    hist(
      outcomes[[col_name]],
      main = paste("Histogramm:", col_name),
      xlab = col_name
    )
  }
}

## shapiro test
shapiro_results <- data.frame(
  outcome = character(),  # Spaltenname
  p_value   = numeric(),    # p-Wert aus dem Test
  stringsAsFactors = FALSE
)

for (col_name in names(outcomes_num)) {
    test_result <- shapiro.test(outcomes_num[[col_name]])
    shapiro_results <- rbind(shapiro_results, data.frame(
        Phenotype = col_name,
        p_value   = test_result$p.value,
        stringsAsFactors = FALSE
      )
    )
}

shapiro_results$normality <- NA

for (i in 1:nrow(shapiro_results)) {
  if (shapiro_results$p_value[i] >= 0.05) {
    shapiro_results$normality[i] <- "normal"
  } else {
    shapiro_results$normality[i] <- "not_normal"
  }
}

table(shapiro_results$normality)

### keine normale datenverteilung: Kruskal Wallis

# ==================================
# Kruskal Wallis
# ==================================

kruskal.result <- kruskal.test(blood_hdl ~ factor(cluster), data=comb_outcomes_cluster10)
kruskal.result


kruskal_results_table <- data.frame(
  outcome = character(), 
  p_value = numeric(), 
  significance = character(), 
  stringsAsFactors = FALSE
)

for(vars in names(comb_outcomes_cluster10[1:23])){
  test_data <- comb_outcomes_cluster10[!is.na(comb_outcomes_cluster10[[vars]]), ] ## spalte an der ich gerade interessiert bin NA's rauslöschen
  if(nrow(test_data) > 0) {
    fml <- as.formula(paste(vars, "~ factor(cluster)"))
    kruskal.result <- kruskal.test(fml, data = test_data)
    if (kruskal.result$p.value < 0.05) {
      sig <- "significant"
    } else {
      sig <- "not significant"
    }
    kruskal_results_table <- rbind(kruskal_results_table, data.frame(phenotype=vars, p_value=kruskal.result$p.value, significance = sig, stringsAsFactors = FALSE))
  }
}


###### DUNN TEST UND BOXPLOTS#######
library(dunn.test)

pdf("Significant_Phenotypes_Boxplots_kmeans10_numerical.pdf", width = 8, height = 8)

for(i in 1:nrow(kruskal_results_table)){
  if(kruskal_results_table$significance[i] != "significant"){
    next #zeilen die nicht signifikant sind werden übersprungen
  }
  pheno <- kruskal_results_table$phenotype[i]
  
  rel_col <- comb_outcomes_cluster10 %>%
    filter(!is.na(.data[[pheno]]))  #### aus relavanter Spalte alle NA rausfiltern (bzw. nur werte lassen die nicht NA sind)
  
  dunn.res <- dunn.test(rel_col[[pheno]], rel_col$cluster, method = "bonferroni")
  
  #save data into dataframe and separate cluster notation -> only keep (filter) data that has p_adjusted value < 0.05
  dunn_df <- data.frame(Comparison = dunn.res$comparisons, p_adjusted = dunn.res$P.adjusted, stringsAsFactors = FALSE ) %>%
    separate(Comparison, into = c("clusterA", "clusterB"), sep = " - ") %>%
    mutate(clusterA = factor(clusterA), clusterB = factor(clusterB)) %>%
    filter(p_adjusted < 0.05)

  if (length(rownames(dunn_df)) == 0){
    next
  }
  
  comparisons_list <- lapply(1:nrow(dunn_df), function(i) {
    c(as.character(dunn_df$clusterA[i]), as.character(dunn_df$clusterB[i]))
  })
  
  annotations_list <- sapply(dunn_df$p_adjusted, function(p) {
    if(p < 0.001) {
      return("***")
    } else if(p < 0.01) {
      return("**")
    } else if(p < 0.05) {
      return("*")
    } else {
      return("ns")
    }
  })
  
  #clustergrößen berechnen für Beschriftung
  cluster_counts <- test_data %>%
    group_by(cluster) %>%
    summarise(n = n())
  
  # my own labels
  new_labels <- cluster_counts %>%
    arrange(cluster) %>%
    mutate(label = paste0(cluster, "\n(n=", n, ")")) %>%
    pull(label)
  
   p <- ggplot(test_data, aes(x = factor(cluster), y = .data[[pheno]])) + ###kommentar
    geom_boxplot() +
    geom_signif(
      comparisons = comparisons_list,
      annotations = annotations_list,
      step_increase = 0.1, 
      map_signif_level = FALSE
    ) +
    scale_y_log10() + ### log10 angewendet um bei 2 ausreißern zB daten besser dazustellen und erkennen
    scale_x_discrete(labels = new_labels) + ### verwendung eigener label
    annotate("text", x = Inf, y = Inf, label = "*: p < 0.05\n**: p < 0.01\n***: p < 0.001",
             hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
    labs(x = "Cluster", y = pheno, title = paste(pheno, " Werte nach Cluster - signifikant nach Dunn Test"))
  
  print(p)
}

dev.off()

#### einzeln ####

test_data <- comb_outcomes_cluster10[!is.na(comb_outcomes_cluster10[["blood_glucose"]]),]

# Führe den Dunn-Test durch, z. B. mit Bonferroni-Korrektur:
dunn.res <- dunn.test(test_data[["blood_hdl"]], test_data$cluster, method = "bonferroni")
print(dunn.res)


#save data into dataframe and separate cluster notation -> only keep (filter) data that has p_adjusted value < 0.05
dunn_df <- data.frame(Comparison = dunn.res$comparisons, p_adjusted = dunn.res$P.adjusted, stringsAsFactors = FALSE ) %>%
  separate(Comparison, into = c("clusterA", "clusterB"), sep = " - ") %>%
  mutate(clusterA = factor(clusterA), clusterB = factor(clusterB)) %>%
  filter(p_adjusted < 0.05)


comparisons_list <- lapply(1:nrow(dunn_df), function(i) {
  c(as.character(dunn_df$clusterA[i]), as.character(dunn_df$clusterB[i]))
})

annotations_list <- sapply(dunn_df$p_adjusted, function(p) {
  if(p < 0.001) {
    return("***")
  } else if(p < 0.01) {
    return("**")
  } else if(p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
})

# boxplots mit signifikanzmerkierungen = zwei cluster die sich signifikant unterscheiden (kleiner p-value = wahrschienlichekeit dass es zufällig ist = so gerung)
ggplot(comb_outcomes_cluster10, aes(x = factor(cluster), y = blood_glucose)) +
  geom_boxplot() +
  geom_signif(
    comparisons = comparisons_list,
    annotations = annotations_list,
    step_increase = 0.1, 
    map_signif_level = FALSE
  ) +
  annotate("text", x = Inf, y = Inf, label = "*: p < 0.05\n**: p < 0.01\n***: p < 0.001",
           hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
  labs(x = "Cluster", y = "blood_glucose", title = paste("blood_glucose Werte nach Cluster - signifikant nach Dunn Test"), shape = "Anzahl Werte")
 #### ende ####

#### CATEGORIAL
outcomes_cat <- select_if(outcomes, is.character)
cat_comb_outcomes_cluster10 <- outcomes_cat
cat_comb_outcomes_cluster10$cluster <- cluster_kmeans

pdf("Distribution_Plot_kmeans10_categorial.pdf", width = 8, height = 8)
for (pheno_name in names(cat_comb_outcomes_cluster10)){
  if (pheno_name == "cluster") next
  plot <- ggplot(cat_comb_outcomes_cluster10, aes(x = factor(cluster), fill = .data[[pheno_name]])) +
    geom_bar() +
    labs(title = paste("Distribution of", pheno_name, "in cluster"),
         x = "cluster",
         y = "number",
         fill = pheno_name) +
    theme_bw()
  
  print(plot)
}
dev.off()

pdf("Significant_Phenotypes_Boxplots_kmeans10_categorial.pdf", width = 8, height = 8)

cat_results_table <- data.frame(
  outcome = character(), 
  p_value = numeric(), 
  method = character(),
  sig = character(),
  stringsAsFactors = FALSE
)

for(pheno in names(cat_comb_outcomes_cluster10)) {
  if(pheno == "cluster") next
  tab <- table(cat_comb_outcomes_cluster10$cluster, cat_comb_outcomes_cluster10[[pheno]])
  
  expected <- chisq.test(tab)$expected #erwartete zellwerte (wichtig für fishers exact beiw eniger als 5)
  #print(chi_result)
  #print(chi_result$observed)
  
  #### Chi Square or Fisher's Exact (if any expected values below 5)
  if (any(expected < 5)){
    test_result <- fisher.test(tab, simulate.p.value = TRUE)
    test <- "Fisher's Exact Test"
  } else {
    test_result <- chisq.test(tab)
    test <- "Chi-Squared Test"
  }
  
  cat_results_table <- rbind(cat_results_table, data.frame(phenotype = pheno, p_value = round(test_result$p.value, 10), method = test, stringsAsFactors = FALSE))
}

#adjusted p-values (Bonferroni-Korrektur und Benjamini-Hochberg-Korrektur)
cat_results_table$p_adjust_bonferroni <- p.adjust(cat_results_table[[2]], method = "bonferroni")
cat_results_table$p_adjust_bh <- p.adjust(cat_results_table[[2]], method = "BH")

for (i in 1:nrow(cat_results_table)){
  print(cat_results_table$p_adjust_bh[i])
  if (cat_results_table$p_adjust_bh[i] < 0.05) {
    print(cat_results_table$p_adjust_bh[i])
    cat_results_table$sig[i] <- "significant"
    print(cat_results_table$sig[i])
  } else {
    cat_results_table$sig[i] <- "not significant"
  }
}
###### ? wie plots? wie pro cluster schauen? dann bleibt nur eine zeile?

if (test_result$p.value < 0.05) {
    sig <- "significant"
    res <- as.data.frame(chi_result$stdres)
    names(res) <- c("cluster", "outcome", "std_resid")
  
    q <- ggplot(res, aes(x = factor(cluster), y = std_resid, fill = outcome)) +
      geom_col(position = "dodge") +  # Balken nebeneinander platzieren
      geom_hline(yintercept = 0, linetype = "dashed") +  # Null-Linie als Referenz+
      #geom_hline(yintercept = cbind(2,-2))
      labs( x = "Cluster",
            y = "Standardized Residuals",
            title = "Chi-Square Standardized Residuals by Cluster",
            fill = names(cat_comb_outcomes_cluster10)[i]) +
      theme_minimal()
    print(q)
  } else {
    sig <- "not significant"
  }
  

dev.off()
  

# Erstelle eine Kreuztabelle
tbl <- table(cat_comb_outcomes_cluster10$cluster, cat_comb_outcomes_cluster10$blood_fasting)
print(tbl)

# Führe den Chi-Quadrat-Test durch:
chi_result <- chisq.test(tbl)
chi_result$expected ##erwartete Häufigkeiten
print(chi_result)

# Fishers Exact Test
fisher_result <- fisher.test(tbl, simulate.p.value = TRUE, B = 1e5)
print(fisher_result)

std_resid <- data.frame(chi_result$stdres)
names(std_resid) <- c("cluster", "outcome", "std_resid")
print(std_resid) ## -2 bis 2 wie erwartet, drunter = deutet auf signifikante unterpräsentation hin - niedriger als erwertet. bei höher = mehr als erwartete häufigkeit

#print(chi_result)

library(ggplot2)

ggplot(std_resid, aes(x = factor(cluster), y = std_resid, fill = outcome)) +
  geom_col(position = "dodge") +  # Balken nebeneinander platzieren
  geom_hline(yintercept = 0, linetype = "dashed") +  # Null-Linie als Referenz
  labs( x = "Cluster",
        y = "Standardized Residuals",
        title = "Chi-Square Standardized Residuals by Cluster",
        fill = "Phenotype") +
  theme_minimal()


ggplot(std_resid, aes(x = factor(cluster), y = std_resid, fill = outcome)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(x = "Cluster", y = "Standardisierte Residuen", title = "Verteilung der standardisierten Residuen pro Cluster") +
  theme_minimal()

ggplot(std_resid, aes(x = factor(cluster), y = std_resid)) +
  geom_boxplot() +
  facet_wrap(~ outcome) +
  labs(x = "Cluster", y = "Frequency", title = "Boxplots der Frequencies pro Cluster für jeden Phenotyp")

# ==================================
# Visualization
# ==================================

data_scaled <- readRDS("/Users/denise/Desktop/BA/data_scaled_and_adjusted.rds")
##### T-SNE #####
# visualizationcluster# visualization mit t-sne
tsne_out <- Rtsne(as.matrix(data_scaled), perplexity = 30, check_duplicates = FALSE, verbose = TRUE)

# In tsne_out$Y finden sich die 2D-Koordinaten
plot(tsne_out$Y, col = cluster_kmeans, pch = 19,
     xlab = "t-SNE 1", ylab = "t-SNE 2",
     main = "Visualisierung der kmm cluster via t-SNE"
)


##### PCA #####
#cluster visualizaiton mit PCA KMEANS https://www.datanovia.com/en/blog/k-means-clustering-visualization-in-r-step-by-step-guide/
cluster_sizes <- table(kmeans_results10$cluster)
kmeans_results10$cluster <- factor(kmeans_results10$cluster, labels = new_labels)
plot <- fviz_cluster(kmeans_results10, data = data_scaled,
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


# gmm
#fviz_cluster(
#  object = list(data = data_scaled, cluster = cluster_assignments),
#  geom = "point",
#  ellipse.type = "convex", 
#  ggtheme = theme_bw()
#)
