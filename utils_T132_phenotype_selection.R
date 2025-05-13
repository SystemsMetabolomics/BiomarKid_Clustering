source("useful_functions.R")
source("utils_phenotype_selection.R")

create_selected_phenotypes_for_later_timepoint <- function(T66, T132, return_overlap_dataset = TRUE){
  T132 <- select_phenotypes(T132)
  T132 <- samples_as_rownames(T132)
  
  ## Refactoring values that are char 
  T132 <- T132 %>%
    mutate(across(where(is.character), as.factor))
  
  T132 <- rownames_as_column(T132)
  
  T66$Sample <- rmv_tp(T66$Sample)
  T132$Sample <- rmv_tp(T132$Sample)
  
  freq <- T66$Sample %in% T132$Sample
  perc <- (sum(freq) / length(T66$Sample)) * 100
  
  T132_valid <- T132[T132$Sample %in% T66$Sample, ] #Nur Samples die in T66 vorkommen
  names_col_na <- names(T132_valid[, colMeans(is.na(T132_valid)) > 0.5])
  
  T132_valid_na <- T132_valid[, colMeans(is.na(T132_valid)) < 0.5]
  
  saveRDS(T132_valid_na, "data/processed/selected_phenotypes_T132_as_in_T66.rds")
  if (return_overlap_dataset){
    return (perc)
  }
}
