library(readxl)
library(writexl)
library(dplyr)
library(magrittr)
library(tidyr)
library(patchwork) # needed to combine plots
library(factoextra) # for PCA (screeplot)
library(ggplot2)
library(openxlsx)

# ==================================
# Read in Data 
# ==================================
# Choose if raw data or clean data (outliers removed)

##### raw data #####
#data <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_NOTwinsorized.xlsx"))
data <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_winsorized.xlsx")) 
selected_phenotypes <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_winsorized.xlsx", sheet = "Selected Phenotype Data"))
all_phenotypes <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_winsorized.xlsx", sheet = "Phenotype Data"))

##### Clean Data #####
data <- as.data.frame(read_excel("/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx")) 
selected_phenotypes <- as.data.frame(read_excel("/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx", sheet = "Selected Phenotype Data"))
all_phenotypes <- as.data.frame(read_excel("/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx", sheet = "Phenotype Data"))

# ==================================
# Data Preparation/Normalization
# ==================================
## set rownames
rownames(data) <- data[,1]
data <- data[,-1]

#str(selected_phenotypes) # to see if they are all num or factorial

## refactoring values that are char
selected_phenotypes <- selected_phenotypes %>%
  mutate(across(where(is.character), as.factor))

## See Data distribution (skewed or not) per metabolite and save as pdf
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
    
    # density plot
    #p_dens <- ggplot(data, aes_string(x = metabolite)) +
      #geom_density(color = border_color, alpha = 0.7) +  # Dichte-Linie
      #ggtitle(paste("Density Plot for", metabolite)) +
      #theme_minimal()
    
    #print(p_dens)
  }
  if (!is.null(pdf_file)) {
    dev.off()
  }
}

#plot_histograms(data, pdf_file= "Metabolites_Raw_Histogramm_wind_Dichtedigramm.pdf")

# -----------------------------
# Log2 Transformation
# -----------------------------
data_log <- log2(data + 1)
#plot_histograms(data_log, pdf_file= "Metabolites_Log_Histogramm_wind_Dichtdiagramm.pdf")

# -----------------------------
# Scaling (Z-Scores)
# -----------------------------
data_scaled <- as.data.frame(scale(data_log))

# ==================================
# Control Data for "Country"
# ==================================
countries <- selected_phenotypes[,1:2]
countries$Sample <- trimws(countries$Sample) #remove spaces

#check if ID's are identical
identical(countries$Sample, rownames(data_scaled)) 

#regression model for each metabolite (each column)
adjusted_metabolites <- data_scaled
for (met in colnames(adjusted_metabolites)) {
  model <- lm(adjusted_metabolites[[met]] ~ countries$Country)
  adjusted_metabolites[[met]] <- residuals(model) 
}

## single metabolite
#model <- lm(adjusted_metabolites$PC.aa.C40.6 ~ countries$Country)
#summary(model)

data_scaled <- adjusted_metabolites

### save preprocessed RDS data for clustering
#saveRDS(data_scaled, file = "data_scaled_and_adjusted.rds")

# ==================================
# PCA to find outliers
# ==================================
pca <- prcomp(data_scaled, scale. = F) # scale = false because data is already scaled
scores <- pca$x  # PCA-transformed data

## scree plot
summary(pca)
fviz_eig(pca, addlabels = TRUE) # nur 30% in ersten 2 Components

# BiPlot (https://www.datacamp.com/de/tutorial/pca-analysis-r)
fviz_pca_var(pca, col.var = "black") # metabolites/Pfeile
fviz_pca_ind(pca, col.ind = "black", labels = T) # nur Punkte (samples)
fviz_pca_biplot(pca, repel = T, col.var = "steelblue", col.ind= "orange") #blau = metabolite, orange = Proben

#pdf("pca_mahalanobis_plot_wind_ca_0.99999.pdf", width = 8, height = 6)

# Visualize first two components
plot(scores[, 1], scores[, 2], main = "PCA", xlab = "PC1", ylab = "PC2", pch = 16)

# Compute Mahalanobis distance (across all Dimensions)
dist_mahal <- mahalanobis(scores, colMeans(scores), cov(scores))
cutoff <- qchisq(0.9999, df = ncol(scores))  # 99% confidence, bei 0.9999999 sieht man nur die 2 eindeitigen outlier vom hc. ## heißt 0.01% liegen nur drüber 
outliers <- which(dist_mahal > cutoff)

# Mark outliers
points(scores[outliers, 1], scores[outliers, 2], col = "red", pch = 19)
text(scores[outliers, 1], scores[outliers, 2], labels = rownames(data_scaled)[outliers],col = "red", cex = 0.8, pos = 3)  # pos = 3: oberhalb des Punktes
print(outliers)  # Outlier indices

#dev.off()


#### REMOVE OUTLIERS & erstelle neue Excel Datei ####
data_rem_outl <- data[-outliers, ]
outcome_rem_outl <- selected_phenotypes[!rownames(selected_phenotypes) %in% outliers, ]
all_pheno_rem_outl <- all_phenotypes[!rownames(all_phenotypes) %in% outliers,]

wb <- createWorkbook()

# create worksheets
addWorksheet(wb, "Metabolomics Data")
addWorksheet(wb, "Selected Phenotype Data")
addWorksheet(wb, "Phenotype Data")

# write data in worksheets
writeData(wb, "Metabolomics Data", data_rem_outl)
writeData(wb, "Selected Phenotype Data", outcome_rem_outl)
writeData(wb, "Phenotype Data", all_pheno_rem_outl)
saveWorkbook(wb, "CHOP_T66_winsorized_no_outliers.xlsx", overwrite = TRUE)