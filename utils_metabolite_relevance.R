source("00_init.R")

### sort rows (nach cluster)
sort_cluster <- function(data, cluster_assignments, k){
  if (identical(rownames(data), rownames(cluster_assignments))){
    sorted_clipped <- data
    sorted_clipped$cluster <- cluster_assignments[[paste0("k_",k)]]
    sorted_clipped <- sorted_clipped[order(sorted_clipped$cluster), ]
  } else {
    stop("Error in sort_cluster: row names of data and cluster assignments are not identical.")
  }
  return (sorted_clipped)}

create_cluster_means <- function(data, cluster_assignments, k){
  data_tmp <- data_scaled
  data_tmp$cluster <- cluster_assignments[[paste0("k_", k)]]
  cluster_means <- aggregate(. ~ cluster,  # everything split into groups according to $cluster
                             data = data_tmp,
                             mean) # mean is calculated for each "group" for all columns
  
  cluster_means <- cluster_means[, !(names(cluster_means) == "cluster")]
  return (cluster_means)
}

rank_metabolites <- function(cluster_means, top_n = 15, bottom_n = 15, filename){
  output <- data.frame(cluster = character(),
                       RankType = character(),
                       Metabolite = character(),
                       Value = numeric()
  )
  
  cluster_means_scaled <- as.data.frame(scale(cluster_means))
  for (cluster in rownames(cluster_means_scaled)){
    
    # create vector of row and its values
    values <- unlist(cluster_means_scaled[cluster,])
    
    top <- sort(values, decreasing = TRUE)[1:top_n]
    top_df <- data.frame(
      Cluster = cluster,
      RankType = "Top",
      Metabolite = names(top),
      Value = as.numeric(top)
    )
    
    bottom <- sort(values, decreasing = FALSE)[1:bottom_n]
    bottom_df <- data.frame(
      Cluster = cluster,
      RankType = "Bottom",
      Metabolite = names(bottom),
      Value = as.numeric(bottom)
    )
    
    output <- rbind(output, top_df, bottom_df)
  }
  
  write.csv(output, file = filename, row.names = FALSE)
}
