### Packages installation for the MILES Reproducibility study ###
## CRAN Packages
# Need "vegan", "dplyr", "tidyr", "ggplot2", "here"

all_cran_pkgs <- c("knitr", "rmarkdown", "BiocManager", "utils", "rvest",
                   "janitor", "here", "tidyverse", "jsonlite", "data.table", 
                   "duckdb", "CDMConnector", "renv", "arrow", "patchwork", 
                   "devtools", "remotes")
install.packages(all_cran_pkgs)

## Phyloseq package
install.packages("BiocManager")
BiocManager::install("phyloseq")

## GitHub Packages - Microshades
remotes::install_github("KarstensLab/microshades", dependencies = TRUE)


