# ==================================
# Create Directories 
# ==================================
directories <- c(
  "data/raw",
  "data/processed",
  "output/plots"
)

invisible(lapply(directories, function(d){
  if (!dir.exists(d)){
    dir.create(d, recursive = TRUE)
    message ("Created directory: ", d)
  }
}))

# ==================================
# Install Packages if required
# ==================================
packages <- c(
  "dplyr",
  "factoextra", # for scree plot PCS
  "ggplot2",
  "magrittr",
  "patchwork", # needed to combine ggplots
  "readxl",
  "tidyr",
  "writexl",
  "aweSOM",
  "cluster",
  "clusterSim",
  "clValid",
  "fpc",
  "heatmaply",
  "kohonen",
  "mclust",
  "NbClust",
  "openxlsx",
  "pheatmap",
  "rJava",
  "subspace",
  "ggalluvial",
  "BiocManager",
  "ComplexHeatmap",
  "RColorBrewer",
  "ggsignif",
  "Rtsne",
  "tibble",
  "ClusterR",
  "devtools",
  "Hmisc",
  "lsr",
  "polycor",
  "circlize"
)


to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

if(length(to_install) > 0){
  install.packages(to_install)
}

# packages that are not available via CRAN
if(!requireNamespace("sgi", quietly = TRUE)){
  devtools::install_github(repo="krumsieklab/sgi", subdir="sgi")
}

if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)){
  BiocManager::install("ConsensusClusterPlus")
}

all_packages <- c(packages, "ConsensusClusterPlus", "sgi")

# Load all packages and remove variables
invisible(lapply(all_packages, library, character.only = TRUE))

rm(all_packages, packages, to_install, directories)
