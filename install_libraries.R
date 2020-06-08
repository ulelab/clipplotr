# Script to install necessary libraries for CLIPplotR
# A. M. Chakrabarti

packages <- c("optparse", "BiocManager", "ggplot2", "ggthemes", "cowplot", "smoother", "zoo", "dplyr", "patchwork")
for(package in packages) {

    if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    message("Installing ", package)
    install.packages(package, character.only = TRUE, repos = "https://cloud.r-project.org")
    
  }

}

biocpackages <- c("rtracklayer", "GenomicFeatures", "ggbio")
for(package in biocpackages) {
  
  if(!suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE))) {
    
    message("Installing ", package)
    BiocManager::install(package)
    
  }  
  
}