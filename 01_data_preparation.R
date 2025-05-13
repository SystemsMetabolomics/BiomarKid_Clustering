source("00_init.R")
source("02_outlier_detection.R")
source("useful_functions.R")
source("utils_plotting.R")
source("utils_phenotype_selection.R")
source("utils_data_preparation.R")

## should outliers be detected and removed?
dummy_set <- FALSE

run_outlier_detection <- FALSE ## no outlier for dummy data
save_exploratory_plots <- FALSE

# ==================================
# Read in Data 
# ==================================
raw_data <- readRDS("dummy_data.rds") 
all_phenotypes <- readRDS("dummy_phenotypes_T66.rds")
phenotypes_T132 <- readRDS("dummy_phenotypes_T132.rds")

## for real dataset (excel sheets)
#raw_data <- as.data.frame(read_excel("CHOP_T66_Denise_winsorized.xlsx")) 
#all_phenotypes <- as.data.frame(read_excel("CHOP_T66_Denise_winsorized.xlsx", sheet = "Phenotype Data"))
#phenotypes_T132 <- as.data.frame(read_excel("CHOP_T132_Denise_winsorized.xlsx", sheet = "Phenotype Data"))


if(run_outlier_detection){
  data_no_outlier <- outlier_detection(raw_data, all_phenotypes, save_outlier_plot = TRUE, create_exp_plots = TRUE)
  data <- data_no_outlier[[1]]
  all_phenotypes <- data_no_outlier[[2]]
  rm(data_no_outlier)
} else {
  data <- raw_data
}

# ==================================
# Phenotype Selection
# ==================================
## Check structure of phenotypes
#str(all_phenotypes) # to see if they are all num or factorial

## Refactoring values that are char 
all_phenotypes <- all_phenotypes %>%
  mutate(across(where(is.character), as.factor))

## select phenotypes (based on visually inspecting correlation matrix see exploratory plots)
selected_phenotypes <- select_phenotypes(all_phenotypes)

# ==================================
# Data Preprocessing
# ==================================
data_log_z <- data_preprocessing(data)

# ==================================
# Adjust Data for Country
# ==================================
clean_data <- adjust_for_country(data_log_z, selected_phenotypes)

# ==================================
# Export RDS
# ==================================
if (run_outlier_detection){
  ### save preprocessed RDS data for clustering
  saveRDS(selected_phenotypes, "data/processed/selected_phenotypes_no_outliers.rds")
  saveRDS(clean_data, file = "data/processed/metabolites_clean_no_outliers.rds")
  saveRDS(phenotypes_T132, file = "data/raw/phenotypes_T132.rds")
} else {
  if(dummy_set){
    saveRDS(data, "data/raw/metabolites_no_outliers.rds")
    saveRDS(all_phenotypes,"data/raw/phenotypes_no_outliers.rds")
    saveRDS(phenotypes_T132, file = "data/raw/phenotypes_T132.rds")
    
    ### save preprocessed RDS data for clustering
    saveRDS(selected_phenotypes, "data/processed/selected_phenotypes_no_outliers.rds")
    saveRDS(clean_data, file = "data/processed/metabolites_clean_no_outliers.rds")
  } else {
    ## raw data
    saveRDS(data, "data/raw/all_metabolites_raw.rds")
    saveRDS(all_phenotypes, "data/raw/all_phenotypes_untouched.rds")
    saveRDS(phenotypes_T132, file = "data/raw/phenotypes_T132.rds")
    
    ### save preprocessed RDS data for clustering
    saveRDS(selected_phenotypes, "data/processed/selected_phenotypes.rds")
    saveRDS(clean_data, file = "data/processed/metabolites_clean.rds")
  }
}

# ==================================
# Exploratory Plots
# ==================================
if (save_exploratory_plots){
  # metabolite histogram raw data
  plot_histograms(raw_data, pdf_file= "output/plots/metabolites_raw_histograms.pdf")
  
  # metbolite histograms log2 transformed data (no outliers if action is chosen)
  plot_histograms(log2(data[,-1]), pdf_file= "output/plots/metabolites_log_histograms.pdf")
  
  #### Correlation plot of numerical phenotypes - Samples as rownames ok because they get removed in function
  create_correlation_plot(all_phenotypes, pdf_file = "output/plots/phenotype_correlation_matrix.pdf")
}


