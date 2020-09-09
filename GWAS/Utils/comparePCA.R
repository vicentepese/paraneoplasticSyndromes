# Import libraries
library(xlsx)
library(readxl)
library(jsonlite)
library(dplyr)
library(stringr)

########### IMPORT #########
setwd("~/Documents/anti-HU")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

# Import Adi PCA
adi.PCA <- read.table('~/Downloads/AutoPCAMatched_LGL1.csv', sep=',', header = TRUE)

# Import list of patients 
vic.PCA <- read.table(settings$file$patList)

# Common IDS
IDS.comm <- merge(adi.PCA, vic.PCA, by.x = 'ID1', by.y = 'V1')
IDs <- IDS.comm$ID1
Dx.comm <- IDS.comm$DX
table(Dx.comm)

# Manual computation
adi.ids <- adi.PCA$ID1
manual.IDs <- c()
for (ID in vic.PCA$V1){
  if (ID %in% adi.ids){
    manual.IDs <- c(manual.IDs, ID)
  }  
}
table(manual.Dx)

# Check IBD filter with aditya 
genome <- read.table('Resources/HU.genome', sep = '\t', header = TRUE)
highIBD.genome <- genome[which(genome$PropIBD > 0.2),]
ID.IBD <- merge(adi.PCA, highIBD.genome, by.x = 'ID1', by.y = 'FID1')
