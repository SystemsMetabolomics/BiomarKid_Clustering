This Bachelor Thesis consists of different scripts. For running this code,
a dummy dataset is included since the real dataset is not publicly distributable.
The dummy set is just random data that follows the structure of the original
data but not data this code was written for. 

Please run 00.init.R first.
This will load all packages required and create the directories needed to save
tables and plots during the analyses. This will create the following structure:

├─ data/                                  # Directory is created by 00_init.R 
│   ├─ raw/                               # Directory is created by 00_init.R 
│   └─ processed/                         # Directory is created by 00_init.R   
├─ 00_init.R               
├─ 01_data_preparation.R
├─ 02_clustering_and_internal_evaluation.R
├─ 03_external_evaluation.R
├─ 04_phenotype_association.R
├─ 05_metabolite_relevance.R
├─ useful_functions.R
├─ utils_data_preparation.R
├─ utils_outlier_detection.R 
├─ utils_phenotype_selection.R
├─ utils_T132_phenotype_selection.R
├─ utils_phenotype_association.R
├─ utils_metabolite_relevance.R
├─ utils_plotting.R
├─ method_exploration.R                 
├─ dummy_data.rds                       # Dummy Data for analyses
├─ dummy_phenotypes_T66.rds             # Dummy Data for analyses
├─ dummy_phenotypes_T132.rds            # Dummy Data for analyses
├─ output/                              # Directory is created by 00_init.R
│   └─ plots/                           # Directory is created by 00_init.R
├─ Supplementary_files/                 # Supplementary Files for Thesis
└─ README.md

All main scripts with XX_script must be run in order, since some of them create 
data needed for the next script. utils-scripts consist of functions and are 
called with source() in the according main scripts. All scripts can only be
run on a specific dataset and do not function as a pipeline for data of 
different structures, since specific phenotypes etc. are selected. 

method_exploration can not be run as a whole, but is rather a file where
different parameters and algorithms must be tested separately. This does not 
create any data needed for the main scripts, but helps chosing the parameters
for the clustering by visually expecting the outcomes of these exploratory 
methods.

01_data_preparation.R -> reads file for analyses in and preprocesses it (log2,
scaling, outlier detection, phenotype selection, adjust data for country)
Metablite input has to be a table of metabolites (columns) and samples (rows) 
or phenotypes (columns) and samples (rows) for phenotype association.
Use dummy_set <- TRUE for dummy set and leave outlier detection on true ->
if not - paths to files need to changed in other scripts.

02_clustering_and_internal_evaluation.R -> here clustering is done on selected
parameters. This results in cluster assignment result tables for k = 2:16
and plots are created for each k and evaluation metric. bootstrapping can be
turned to TRUE (but takes a while to fully run)
This file script automatically reads in file created from 01.

03_external_evaluation.R -> creates evaluation plot for all methods to compare.
creates ARI plot and alliuvial plots for visuall inspection.

04_phenotype_association.R -> timepoint for which phenotype association should be 
done has to be chosen here (T66 or T132). Also select k_values association should
be done on (between 2:16) for association tables and a k (must be in 
k_values) to create plots for.

For SGI analysis: 05 SGI -> TRUE (only works for T66)

05_metabolite_relevance.R -> creates heatmaply.html, a circle mean heatmap and 
ranks metabolites that drive clusters most (top10 and bottom10). Choose k here. 



