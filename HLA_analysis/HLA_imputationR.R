# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(parallel)
library(HIBAG)

## DESCRIPTION ##
# This function takes a PLINK format master file (including both cases and controls)
# And computes HLA imputation. HLA calls are stored individually for each locus

########## IMPORT ##########
setwd("~/Documents/paraneoplasticSyndromes")

# Import settings
settings <- jsonlite::fromJSON('HLA_analysis/settingsHLA.json')

######### COLLAPSE PRE-TAINED MODELS #########

# Load main model 
mainModel <- get(load(settings$file$PMRA_HLA_model))

# LOAD DRB3-5 models and remove sample.id
modelDRB3 <- get(load(settings$file$PMRA_DRB3_model)); modelDRB3[['sample.id']] <- NULL
modelDRB4 <- get(load(settings$file$PMRA_DRB4_model)); modelDRB4[['sample.id']] <- NULL
modelDRB5 <- get(load(settings$file$PMRA_DRB5_model)); modelDRB5[['sample.id']] <- NULL

# Add DRB3-5 to main file 
mainModel[['DRB3']] <- modelDRB3; mainModel[['DRB4']] <- modelDRB4; mainModel[['DRB5']] <- modelDRB5;

# Re-save 
save(mainModel, file = settings$file$PMRA_HLA_model)


########## HLA IMPUTATION ############

# List of models 
models <- list.files(settings$directory$HLA_trained_models)

# Load pre-fit model and comvert to hlaMODEL
model.list <- get(load(settings$file$PMRA_HLA_model))
hla.id <- names(model.list)

# Import file
gname <- settings$file$HU
yourgeno <- hlaBED2Geno(bed.fn=paste(gname, ".bed", sep = ''), fam.fn=paste(gname, ".fam", sep='')
                        , bim.fn=paste(gname, ".bim", sep=''), assembly = 'hg19')
summary(yourgeno)

# Make cluster 
cl <- makeCluster(10)

# Make predictions
for (locus in hla.id){
  model.hla <- hlaModelFromObj(model.list[[locus]])
  summary(model.hla)
  pred.guess <- predict(model.hla, yourgeno, type="response+prob", nclassifiers=100, cl=cl)
  save(pred.guess, file = paste(settings$directory$HLA_predictions, paste('HLA', locus, sep = ''), '.RData', sep= ''))
}


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
write.table(data, file = settings$file$HLA_total, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ',')


