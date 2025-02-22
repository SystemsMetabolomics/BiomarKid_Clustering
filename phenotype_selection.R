library(readxl)
library(dplyr)
library(corrplot)
library(Hmisc)
library(openxlsx)
library(lsr)
library(polycor)

phenotypes <- as.data.frame(read_excel("/Users/denise/Desktop/CHOP_T66_Denise_winsorized.xlsx", sheet = "Phenotype Data"))
rownames(phenotypes) <- phenotypes[,1]
phenotypes <- phenotypes[,-1]

# Entfernen von Spalten mit mehr als 50% fehlenden Werten
phenotypes <- phenotypes[, colMeans(is.na(phenotypes)) < 0.5] ##includes all columns that are just NA's

# remove columns that all have the same value - weil sie keine Informationen über Zusammenhänge zwischen den outputs zeigen
unique_count <- sapply(phenotypes, function(x) length(unique(x[!is.na(x)]))) 
unique_cols <- c()

for (col in names(unique_count)) {
  if (unique_count[[col]] == 1) {
    unique_cols <- c(unique_cols, col)
  }
}

phenotypes <- phenotypes[, !(colnames(phenotypes) %in% unique_cols)]

## refactoring values that are char
phenotypes <- phenotypes %>%
  mutate(across(where(is.character), as.factor))

## function for corrplot https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
my_corrplot <- function(matrix, plot_title = NULL){
  corrplot(matrix,
           order = "hclust",
           diag = FALSE,        # damit diagonale die ja 1 überall ist nicht angezeigt wird
           type = "upper",      # nur obere, reche hälfte anzeigen (links ist ident)
           tl.cex = 0.6,          # Labels kleiner
           tl.col = "black",     # Labels Farbe
           tl.srt = 45,              # Labels drehen zum besseren lesen
           addCoef.col = "black",      # Farbe der Korrelationskoeffizienten
           number.cex = 0.2,           # Größe der Korrelationskoeffizienten
           title = plot_title
           )
}

### Correlation Matrix for numerical only
phenotypes_num <- select_if(phenotypes, is.numeric)
cor_matrix <- cor(phenotypes_num, method = "pearson", use = "pairwise.complete.obs") # tried spearman auch aber ist im großen und ganzen fast ident

my_corrplot(cor_matrix, "Correlation Matrix numerical phenotypes")

### Correlation Matrix for categorical with cramers V https://search.r-project.org/CRAN/refmans/confintr/html/cramersv.html
phenotypes_cat <- phenotypes[, sapply(phenotypes, function(x) is.factor(x))] # filter out factors (not numerics)

# create empty matrix to fill later
cramers_matrix <- matrix(NA, nrow = length(phenotypes_cat), ncol = length(phenotypes_cat), dimnames = list(colnames(phenotypes_cat), colnames(phenotypes_cat)))

# loop through all pairs of variables
for(i in 1:(length(phenotypes_cat) - 1)) {
  for(j in (i + 1):length(phenotypes_cat)) {
    #contingency table
    tab <- table(phenotypes_cat[[i]], phenotypes_cat[[j]])
    cv <- cramersV(tab)
    #save values in matrix
    cramers_matrix[i, j] <- cv
    cramers_matrix[j, i] <- cv
  }
}

## warnings() !!!!!

diag(cramers_matrix) <- 1 #Diagonale auf 1 setzen
cramers_matrix[is.na(cramers_matrix)] <- 0 #NA's auf 0 setzen

my_corrplot(cramers_matrix, "Correlation Matrix categorial")

###### BOTH with polycor
hetcor_result <- hetcor(phenotypes)
mixed_matrix <- hetcor_result$correlations
mixed_matrix[is.na(mixed_matrix)] <- 0

my_corrplot(mixed_matrix, "Correlation Matrix all phenotypes")

### manually looked at all correlations and the phenotypes, and which one we are interested in or which might be useful.
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

#### selected phenotypes in excel sheet mit neuen selected phenotypes erstezen
selected_phenotypes <- selected_phenotypes[-outliers, ] ## outliers von davor entfernen

wb <- loadWorkbook("/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx")
removeWorksheet(wb, "Selected Phenotype Data") ## altes Sheet entfernen
addWorksheet(wb, "Selected Phenotype Data")
writeData(wb, "Selected Phenotype Data", selected_phenotypes)
saveWorkbook(wb, "/Users/denise/Desktop/BA/CHOP_T66_winsorized_no_outliers.xlsx", overwrite = TRUE)