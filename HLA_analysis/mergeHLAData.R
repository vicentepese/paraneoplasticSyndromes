# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(HIBAG)
library(parallel)
library(corrplot)
library(randomForest)
library(xlsx)

## DESCRIPTION ##
# Merge HLA Data 

########## IMPORT ##########
setwd("~/Documents/paraneoplasticSyndromes")

# Import settings
settings <- jsonlite::fromJSON('HLA_analysis/settingsHLA.json')

########## COLLAPSE #########
# Collapse all HLA locus calls in a single file 

# Get prediction files 
path <- settings$directory$HLA_predictions
pred.files <- list.files(settings$directory$HLA_predictions)

# Load first file and get list of samples
init.data <- get(load(paste(path, pred.files[1], sep = '')))
sample.id <- pred.guess$value$sample.id

# Get GWASIDs and append create datafrane
changeIdFormat = function(GWASID){
  idpre <- strsplit(GWASID, '_')[[1]]
  return(strsplit(idpre[length(idpre)], '\\.')[[1]][1])
}
GWASIDs <- sapply(sample.id, changeIdFormat)
data <- data.frame(sample.id = sample.id, GWASID=GWASIDs, row.names = NULL)

# Collapse all HLA calls in a single dataframe 
for (file in pred.files){
  temp <- get(load(paste(path, file, sep = '')))
  name <- strsplit(file, '\\.')[[1]][1]
  data[paste(name, 'A1', sep = '_')] <- as.factor(temp$value$allele1)
  data[paste(name, 'A2',sep = '_')] <- as.factor(temp$value$allele2)
}

# Write 
write.table(data, file = settings$file$HLA_total, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
