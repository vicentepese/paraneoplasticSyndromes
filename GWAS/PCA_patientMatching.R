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
library(zeallot)

############ DESCRIPTION ############
# Case-control matching using Principal Components including all additional subjects. 
# Euclidian distance between each subject and the other is computed using the first two PCs. 
# The closest number of controls as defined in settingsGWAS are taken 

############ INITIALIZATION ###########

# Set working directory 
setwd("~/Documents/paraneoplasticSyndromes/")

# Import settings 
settings <- jsonlite::fromJSON('GWAS/settingsGWAS.json')

# Load PCs
PCs <- read.table(settings$file$PCs_all_plates, header = TRUE)

# Load list of cases
cases <- read.table(settings$file$GWASIDsCases) %>% unlist()

############ PREPROCESS PCS ##############
# Add diagnosis to dataset 

# Convert to GWASID 
convert2GWASID=function(id){
  if (grepl('Stanford', id)){
    return(id %>% strsplit('\\.') %>% unlist() %>% .[1] %>% strsplit('_') %>% unlist() %>% .[4])
  }else{
    return(id)
  }
}
PCs$GWASID <- sapply(PCs$FID, convert2GWASID) %>% unlist()

# Get diagnosis
getDx=function(GWASID, cases){
  if (GWASID %in% cases){
    return (1)
  } else{
    return (0)
  }
}
PCs$Dx <- sapply(PCs$GWASID, getDx, cases) %>% unlist()

############# EUCLIDEAN DISTANCE ############

# Separate datasets 
PCs.cases <- PCs %>% filter(Dx == 1)
PCs.controls <- PCs %>% filter(Dx == 0)

# Compute Euclidean distance for each case 
matched.Controls <- c()
for (case in PCs.cases$GWASID){
  
  # Get PCs
  case.data <- PCs.cases %>% filter(GWASID == case)
  PC12.cases <- c(case.data$PC1, case.data$PC2)
  
  PCdist.df <- data.frame()
  for (control in PCs.controls$GWASID){
    
    # Get Pcs 
    control.data <- PCs.controls %>% filter(GWASID == control)
    PC12.controls <- c(control.data$PC1, control.data$PC2)
    
    # Compute euclidean distance 
    PCdist.df <- rbind(PCdist.df,
                    data.frame(GWASID = control.data$GWASID, 
                               PCdist= dist(x = rbind(PC12.cases, PC12.controls), 
                                            method = 'euclidean') %>% .[1]))
  }
  
  # Sort by dist 
  PCdist.df <- PCdist.df[order(PCdist.df$PCdist),]
  
  # Append GWASD ID of matched controls 
  matched.Controls <- c(matched.Controls, c(PCdist.df$GWASID[1:settings$ControlCaseRatio]))
}

# Get unique GWASIDs 
matched.Controls <- unique(matched.Controls)

# Save matched controls 
write.csv(matched.Controls, file = settings$file$matchedControls, quote = TRUE, row.names = FALSE)

##################### PLOT ###################

# Load matched controls 
tst <- read.csv(settings$file$matchedControls) %>% unlist()

# Subset data 
data.plt <- rbind(PCs.cases, PCs.controls %>% filter(GWASID %in% matched.Controls))

#

# Plot 
data.plt %>%
  ggplot(aes(x = PC1, y = PC2, color = as.factor(Dx))) + geom_point()
data.plt %>% filter(Dx == 1) %>% 
  ggplot(aes(x = PC1, y = PC2)) + geom_point()
data.plt %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(Dx))) + geom_point() +xlim(0,0.007) + ylim(-0.005, 0.017)

         