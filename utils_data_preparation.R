# ==================================
# Data Preprocessing
# ==================================
data_preprocessing <- function (data) {
  ### log2 transformation
  data_log <- data
  data_log[, -1] <- log2(data_log[, -1])
  
  ### scaling (z-scores)
  data_log_z <- data_log
  data_log_z[, -1] <- as.data.frame(scale(data_log_z[, -1]))
  
  return (data_log_z)
}

### control data for country
adjust_for_country <- function (data, phenotypes){
  countries <- phenotypes %>%
    select(Sample, Country)
  #countries$Sample <- trimws(countries$Sample) #remove spaces
  
  merged_data <- merge(data, countries, by = "Sample", sort = FALSE)
  merged_data <- samples_as_rownames(merged_data)
  
  #regression model for each metabolite (column)
  adjusted_metabolites_merged <- merged_data
  for (met in colnames(adjusted_metabolites_merged)) {
    if (met != "Country") {
      model <- lm(adjusted_metabolites_merged[[met]] ~ adjusted_metabolites_merged$Country)
      adjusted_metabolites_merged[[met]] <- residuals(model)
    }
  }
  
  ## single metabolite
  #model <- lm(adjusted_metabolites$PC.aa.C40.6 ~ countries$Country)
  #summary(model)
  
  ## clean data without last column (country from merged)
  data_clean <- adjusted_metabolites_merged[, -ncol(adjusted_metabolites_merged)]
  return (data_clean)
}
