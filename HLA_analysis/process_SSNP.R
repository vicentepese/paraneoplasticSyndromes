# Import libraries
library(xlsx)
library(readxl)
library(jsonlite)
library(dplyr)
library(stringr)

########### IMPORT #########
setwd("~/Documents/paraneoplasticSyndromes")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

# Import association file 
assoc.data <- read.table("./Resources/Plate77_121.assoc.logistic", header = FALSE, sep = '')
colnames(assoc.data) <- c( "CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P")
assoc.data <- assoc.data[which(assoc.data$TEST == 'ADD'),]
assoc.data <- assoc.data[order(assoc.data$CHR, assoc.data$BP),]

########## PRE-PROCES ###########

# Get low p-values SNPS
assoc.data['log10P'] <- -log10(assoc.data$P)
assoc.data.lowP <- assoc.data[which(assoc.data$log10P > 5 ),]

# Compute differences between SNPs within each chromosome
CHRs <- unique(assoc.data.lowP$CHR)
dif <- c()
for (chr in CHRs){
  
  # Get data from chromsomes 
  assoc.data.CHR <- assoc.data.lowP[which( assoc.data.lowP$CHR == chr),]
  
  # Compute differences 
  dif <- c(dif, 0)
    if (nrow(assoc.data.CHR) > 1){
    for (i in c(2:nrow(assoc.data.CHR))){
      dif <- c(dif, assoc.data.CHR$BP[i] - assoc.data.CHR$BP[i-1])
    }
  } 
} 

# Append to data
assoc.data.lowP['BPdiff'] <- dif
idxs <- which(assoc.data.lowP$BPdiff < 100000 & assoc.data.lowP$BPdiff > 0)
assoc.data.SameRegion <- list()
for (i in c(1:length(idxs))){
  idx <- idxs [i]
  temp <- assoc.data.lowP[c(idx-1, idx),]
  snp1 <- as.character(temp$SNP[1]); snp2 <- as.character(temp$SNP[2])
  pVal1 <- temp$log10P[1]; pVal2 <- temp$log10P[2]
  bpdist <- temp$BPdiff[2]; chr <- temp$CHR[1]
  assoc.data.SameRegion[[i]] <- c(chr, snp1, snp2, pVal1, pVal2, bpdist) 
}
assoc.data.SameRegion <- as.data.frame(do.call(rbind, assoc.data.SameRegion))
colnames(assoc.data.SameRegion) <- c('CHR', 'SNP1', 'SNP2', 'PVAL1', 'PVAL2', 'DBDIST')
write.xlsx(assoc.data.SameRegion, paste('Outputs/BPdists_', settings$diagnosis, '.xls', sep = ''), col.names = TRUE, row.names = FALSE)
write.xlsx(assoc.data.lowP, paste('Outputs/BPdists_', settings$diagnosis, '_lowP.xls', sep = ''), col.names =  TRUE, row.names = FALSE)
