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
library(readxl)

## DESCRIPTION ##
# This function loads the HLA calls for diseases with low sample sizes
# and creates a file that allows visual exploration of the results 

########## IMPORT ##########
setwd("~/Documents/paraneoplasticSyndromes")

# Import settings
settings <- jsonlite::fromJSON('HLA_analysis/settingsHLA.json')


########## HLA ANALYSIS FUNCTION ##########

HLA_analysis=function(Dx.samples, Dx.GWASIDs, Ds){
  files = list.files(settings$directory$HLA_predictions)
  HLA_calls = list()
  for (f in files){
    
    # Import file 
    load(paste(settings$file$HLA_predictions, f, sep = ''))
    
    # Change ID format to match GWAS IDs
    sample.IDs <- pred.guess$value$sample.id
    changeIdFormat = function(GWASID){
      idpre <- strsplit(GWASID, '_')[[1]]
      return(strsplit(idpre[length(idpre)], '\\.')[[1]][1])
    }
    GWASIDs <- lapply(sample.IDs, changeIdFormat)
    pred.guess$value['GWASIDs'] <- unlist(GWASIDs)
    
    # Get HLA calls from the locus 
    HLA.data <- pred.guess$value[which(pred.guess$value$GWASIDs %in% Dx.GWASIDs),] 
    colnames(HLA.data) <- c('sample.id', paste('HLA_',pred.guess$locus, '_1', sep = ''), 
                            paste('HLA_', pred.guess$locus, '_2', sep = ''),
                            paste('prob_HLA_', pred.guess$locus, sep = ''), 
                            paste('matching_HLA_', pred.guess$locus, sep = ''), 
                            'GWAS_ID')
    
    # Add HLA calls to list
    HLA_calls[colnames(HLA.data)[2:5]] = HLA.data[colnames(HLA.data)[2:5]]
  }
  
  # Convert to dataframe 
  HLA_calls.df <- data.frame(matrix(unlist(HLA_calls), ncol = length(HLA_calls)))
  colnames(HLA_calls.df) <- names(HLA_calls)
  HLA_calls.df <- cbind(data.frame(GWAS.ID = HLA.data$GWAS_ID), HLA_calls.df)
  summary(HLA_calls.df)
  
  # Create Excel file
  write.table(HLA_calls.df, file = paste(settings$directory$HLA_Output_smallSample, Ds, '.xlsx'), 
              sep = ',', quote = FALSE, col.names = TRUE, row.names = FALSE)
}

############ HLA ANALYSES OF SMALL SAMPLE SIZES #############

# Import paraneoplastic syndrome cases data
paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
Dxs <- table(paraneo.cases$Dx)

# Create list of diseases based on diagnosis
Dss <-c('AK5', 'DNER','GABAbr', 'Ma2')
Dxs <- Dxs[Dss]
Dxs

# Compute analysis of small samples
for (Ds in Dss){
  Dx.samples <- paraneo.cases[which(paraneo.cases$Dx == Ds),]
  Dx.GWASIDs <- Dx.samples$`GWAS ID`
  
  # Compute HLA analysis 
  HLA_analysis(Dx.samples, Dx.GWASIDs, Ds)
}
