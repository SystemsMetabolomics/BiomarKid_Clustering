# make sample column rownames
samples_as_rownames <- function(data){
  if (! any(colnames(data) == "Sample")) {
    stop("There is no column 'Sample'")
  }
  
  col <- data[["Sample"]]
  
  if (any(duplicated(col))) {
    stop("Values in 'Sample' have to be unique in order to use them as rownames.")
  }
  
  rownames(data) <- col
  data$Sample <- NULL
  return(data)
}

rownames_as_column <- function(data){
  data$Sample <- rownames(data)
  return(data)
}

# removes timepoint stamp from sample name
rmv_tp <- function(data_col){
  data_col <- sub("_[^_]+$", "", data_col)
  return (data_col)
}
# _ -> findet den letzten Unterstrich.
#[^_]+ -> findet Folge von Zeichen die keine Unterstriche sind (nach letztem Unterstrich)
#$ -> "T24" am Ende des Strings steht.bestätigt dass Ausdruck am Ende steht



