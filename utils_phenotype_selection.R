library(readxl)
library(dplyr)
library(corrplot)
library(Hmisc)
library(openxlsx)
library(lsr)

## introduction for corrplot https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
my_corrplot <- function(matrix){
  corrplot(matrix,
           order = "hclust",
           diag = FALSE,        # damit diagonale die ja 1 überall ist nicht angezeigt wird
           type = "upper",      # nur obere, reche hälfte anzeigen (links ist ident)
           tl.cex = 0.5,          # Labels kleiner
           tl.col = "black",     # Labels Farbe
           tl.srt = 45,              # Labels drehen zum besseren lesen
           #addCoef.col = "black",      # Farbe der Korrelationskoeffizienten
           #number.cex = 0.2,           # Größe der Korrelationskoeffizienten
           title = "Phenotype Correlation Matrix",
           cex.main = 0.5
  )
}


create_correlation_plot <- function(phenotypes, pdf_file = NULL){
  phenotypes <- samples_as_rownames(phenotypes)
  
  # remove phenotypes with more than 50% NA's
  phenotypes <- phenotypes[, colMeans(is.na(phenotypes)) < 0.5] 
  
  # remove phenotypes that all have the same value (because they are constant and don't show any information for correlation)
  phenotypes <- phenotypes %>%
    select(where(~ length(unique(na.omit(.))) > 1))
  
  # create correlation matrix for numerical phenotypes
  phenotypes_numerical <- select_if(phenotypes, is.numeric)
  cor_matrix <- cor(phenotypes_numerical, method = "pearson", use = "pairwise.complete.obs") 
  
  if (!is.null(pdf_file)) {
    pdf(pdf_file, height = 20, width = 5)
  }
  
  # create correlation plot
  my_corrplot(cor_matrix)
  
  if (!is.null(pdf_file)) {
    dev.off()
  }
}

# manual visualization of plot and categorical phenotypes has led to these selected phenotypes
select_phenotypes <- function (phenotypes){
  keep_phenotypes <- c("Sample",
                       "Country", 
                       "season", 
                       "gender", 
                       "formula", 
                       "fm_sl_perc", 
                       "cbmi", 
                       "bmicatwho", 
                       "blood_fasting", 
                       "blood_hdl", 
                       "blood_ldl", 
                       "blood_chol", 
                       "blood_trig", 
                       "blood_urea", 
                       "blood_creatin", 
                       "blood_uricacid", 
                       "blood_insulin", 
                       "blood_proinsulin", 
                       "blood_leptin", 
                       "blood_igf1",
                       "blood_igf_1_free",
                       "blood_hmw_adiponectin",
                       "blood_apo_b",
                       "blood_apo_ai",
                       "blood_hcrp",
                       "homa_index",
                       "tg_hdl_ratio",
                       "blood_glucose",
                       "PA_norm",
                       "MVPA_norm"
  )
  
  selected_phenotypes <- phenotypes[, keep_phenotypes]
  return(selected_phenotypes)
}
