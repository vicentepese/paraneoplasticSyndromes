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

########## IMPORT #########

# Set working directory 
setwd('~/Documents/paraneoplasticSyndromes')

# Import settings 
settings <- jsonlite::fromJSON('HLA_analysis/settingsHLA.json')

# Import HLA data 
data <-read.table(file = settings$file$HLA_total, header = TRUE)

########## PRE-PROCESS DATA ##########

preprocessData=function(cases, controls){
  
  # Import subject list 
  subjects.ID <- read.table(settings$file$patList, colClasses = c('character', 'NULL')) %>% unlist()
  changeIdFormat = function(GWASID){
    idpre <- GWASID %>% strsplit('_')  %>% unlist() 
    return(strsplit(idpre[length(idpre)], '\\.')[[1]][1])
  }
  subjects.GWASID <- sapply(subjects.ID, changeIdFormat, USE.NAMES = FALSE)
  
  # Filter out patients
  data.filt <- data[which(data$sample.id %in% subjects.ID),]
  
  # Add diagnosis to data
  Dx <- c()
  for (subj in as.character(data.filt$GWASID)){
    if (subj %in% controls){
      Dx <- c(Dx, 0)
    } else if (subj %in% cases) {
      Dx <- c(Dx, 1)
    } else{
      Dx <- c(Dx, -9)
    }
  }
  
  # Add to datframe 
  data.filt['Dx'] <- Dx
  
  # Compute cases and controls data 
  data.cases <- data.filt %>% filter(Dx == 1)
  data.controls <- data.filt %>% filter(Dx == 0)
  
  # Initialize
  locus <- settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[1]
  c(allele1, allele2) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
  allele = settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[2]
  
  # Filter based on experiment
  switch(settings$experiment, 
         '1'={
           #  Relative predisposition analysis: considering DR7 positive cases and controls both homo and heterozygous. 
           #    Allele count in both 
           
           # Get DR7 positive cases and controls 
           data.cases <- data.cases %>% .[,c(allele1,allele2)] %>% filter(get(allele1) == allele | get(allele2) == allele)
           
           # Count number of the allele of Interest in data controls
           alleleCount <- data.controls %>% filter(get(allele1) == allele | get(allele2) == allele) %>% nrow() * 2
           
         },
         
         '2'={
           #  Heterozygote comparison: remove allele (e.g., DRB1*07:01) positive from controls, and allele homozygous from cases.
           # Compare the othe allele in cases-hetero with the number of alleles in controls.
           
            # Get allele data in controls and remove it 
           data.controls <- data.controls[,c(allele1, allele2)] %>% filter(get(allele1) != allele & get(allele2) != allele)
           
           # Get heterozygous cases
           carriers <- data.cases %>% filter(get(allele1) == allele | get(allele2) == allele)
           carriers.hetero.length <- carriers[,c(allele1, allele2)] %>% apply(1, unique) %>% lapply(length) %>% unlist()
           carriers <- carriers[which(carriers.hetero.length == 2),]
           data.cases <- data.cases[which(data.cases$GWASID %in% carriers$GWASID),]
           
           # alleleCount = 0 as they're not removed fromt controls 
           alleleCount <- 0 
         })
  
  # Return
  return(list(data.cases, data.controls, data.filt, alleleCount))
  
}

######### ALLELE / CARRIER FREQUENCY #######

computeACFREQ = function(data, locus, Dx){
  # Compute heterozigous, homozigous, or absence count. 
  # Compute allele frequency, count and total.
  # Compute carrier frequency, count and total.
  
  # Allele frequency 
  A1 <- locus %>% paste('_A1', sep = ''); A2 <- locus %>% paste('_A2', sep = '')
  alleles <- list(data[, A1], data[, A2]) %>% unlist() %>% droplevels()
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  alleles.df <- data.frame(allele = alleles.count %>% names(), 
                           alleleCount = alleles.count %>% as.vector(),
                           alleleFreq = alleles.freq %>% as.vector(), 
                           alleleTotal = nrow(data)*2)
  
  
  
  # Heterozigous, homozigous, and absence count
  carriers <- data[,c(A1,A2)]
  A0 <- c(); A1 <- c(); A2 <- c(); 
  for (A in levels(alleles)){
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(alleles), A0 = A0, A1 = A1, A2 = A2)
  
  # Merge 
  ACFREQ.df <- merge(HH.data, alleles.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',6), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',6), sep = ''))
            return(ACFREQ.df)
          }
  )
}
######## CHI2 / FISHER EXACT FUNCTION #########
# Compute Chi2 and Fisher exact function on allele and carrier 
# count between cases and controls 

computeChi2 = function(ACFREQ.df){
  
  # Initialize 
  chi2.pval.alleles <- c() ; chi2.pval.carriers <- c();
  chi2.OR.alleles <- c(); chi2.OR.carriers <- c();
  fishers.pval.alleles <- c(); fishers.pval.carriers <- c()
  fishers.OR.alleles <- c(); fishers.OR.carriers <- c()
  fishers.UI.alleles <- c(); fishers.UI.carriers <- c()
  fishers.LI.alleles <- c(); fishers.LI.carriers <- c()
  
  # For each allele in the locus 
  for (A in ACFREQ.df$allele){
    
    ## Allele frequency 
    # Create contingency table  and run chi2 test
    allele.data <- ACFREQ.df %>% filter(allele == A)
    cont.table.allele <- matrix(c(allele.data$alleleCountCase,
                                  allele.data$alleleCountControl,
                                  allele.data$alleleTotalCase - allele.data$alleleCountCase ,
                                  allele.data$alleleTotalControl - allele.data$alleleCountControl - alleleCount),
                                nrow = 2)
    # Tests
    chi2.allele.res <- chisq.test(x = cont.table.allele)
    fishers.allele.res <- fisher.test(cont.table.allele)
    # Chi2
    chi2.pval.alleles <- c(chi2.pval.alleles, chi2.allele.res$p.value) ;
    chi2.OR.alleles <- c(chi2.OR.alleles, chi2.allele.res$statistic %>% unname())
    # Fishers' exact test 
    fishers.pval.alleles <- c(fishers.pval.alleles, fishers.allele.res$p.value);
    fishers.OR.alleles <- c(fishers.OR.alleles, fishers.allele.res$estimate %>% unname())
    fishers.LI.alleles <- c(fishers.LI.alleles, fishers.allele.res$conf.int[1])
    fishers.UI.alleles <- c(fishers.UI.alleles, fishers.allele.res$conf.int[2])
    
  }
  stats.data <- data.frame(allele = ACFREQ.df$allele, 
                           FishersAllelePVAL = fishers.pval.alleles, FishersAlleleOR = fishers.OR.alleles, 
                           FishersAlleleLI = fishers.LI.alleles, FishersAlleleUI = fishers.UI.alleles,
                           ChiAllelePVAL = chi2.pval.alleles, ChiAlleleOR = chi2.OR.alleles)
  stats.data[is.na(stats.data)] <- 0
  
  return(stats.data)
  
}

######### CHI2 HOMOZYGOSITY / HETEROZYGOSITY HLA ANALYSIS #######

diseases <- c('LGI1')
locus <- paste('HLA', settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[1], sep = '')

for (DS in diseases){
  
  # Get patient information
  system(paste('Rscript --vanilla getPatInfo.R', DS, sep = ' '))
  
  # Import controls and cases, and merge
  cases <- read.table(settings$file$GWASIDsCases) %>% unlist()
  controls <- read.table(settings$file$GWASIDsControls) %>% unlist()
  
  # Pre-preprocess data 
  c(data.cases, data.controls, data.filt, alleleCount) %<-% preprocessData(cases, controls) %>% unlist() %>% invisible()
  
  # Compute allele, carrier, and homo/hetero count in cases and controls 
  ACFREQ.cases <- computeACFREQ(data = data.cases, locus = locus, Dx = 'case' )
  ACFREQ.controls <- computeACFREQ(data = data.controls, locus = locus, Dx = 'control')
  ACFREQ.df <- merge(ACFREQ.cases, ACFREQ.controls, by = "allele")
  
  # Compute Chi2
  stats.data <- computeChi2(ACFREQ.df = ACFREQ.df)
  
  # Create dataframes 
  HLA.alleles.df <- merge(stats.data[,c(1,which(grepl('Allele', colnames(stats.data))))],
                          ACFREQ.df[,c(which(grepl(paste(c('allele'), collapse = '|'), colnames(ACFREQ.df))))],
                          by = 'allele')
  
  # Write to excel output
  switch (settings$experiment,
          "1" ={
            write.xlsx(x = HLA.alleles.df, file = paste(settings$directory$HLA_Output_HomoHetero,
                                                        '/HLA_Relative_Predisposition_', DS, '_', settings$alleleOfInterest,'.xlsx', sep = ''),
                       col.names = TRUE, row.names = FALSE, append = TRUE)
          },
          "2" = {
            write.xlsx(x = HLA.alleles.df, file = paste(settings$directory$HLA_Output_HomoHetero,
                                                        '/HLA_Heterozygote_Comparison_', DS, '_', settings$alleleOfInterest,'.xlsx', sep = ''),
                       col.names = TRUE, row.names = FALSE, append = TRUE)
          })

  
}
