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
if(settings$allPlates == 1){
  patientList <- read.table(settings$file$patList_allplates) %>% unlist()
  data <-read.table(file = settings$file$HLA_calls_all_plates, header = TRUE, sep = ',')
  data <- data %>% filter(GWASID %in% patientList)
}else{
  data <-read.table(file = settings$file$HLA_total, header = TRUE, sep = ',')
}


########## PRE-PROCESS DATA ##########

preprocessData=function(data, settings){
  
  if (settings$allPlates != 1){
    
    # Import controls and cases, and merge
    cases <- read.table(settings$file$GWASIDsCases) %>% unlist()
    controls <- read.table(settings$file$GWASIDsControls) %>% unlist()
    
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
    
  } else{
    data.filt <- data
  }

  # Compute cases and controls data 
  data.cases <- data.filt %>% filter(Dx == 1)
  data.controls <- data.filt %>% filter(Dx == 0)
  
  # Initialize
  locus <- settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[1]
  c(allele1, allele2) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
  allele = settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[2]
  
  #  All controls, heterozygous cases
  # Get heterozygous positive cases
  data.cases <- data.cases %>% .[,c(allele1,allele2)] %>% filter((get(allele1) == allele | get(allele2) == allele) & 
                                                                   !(get(allele1) == allele & get(allele2) == allele))
  # Remove other allele 
  alleleCount <- 0
  for (A in settings$allele2Remove){
    # Locus and allele
    locusREM <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    c(allele1Rem, allele2Rem) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
    alleleRem = settings$allele2Remove %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Count allele in controls 
    alleleCount <- alleleCount + data.controls[,allele1Rem] %>% table() %>% .[alleleRem] %>% unname() +
      data.controls[,allele2Rem] %>% table() %>% .[alleleRem] %>% unname()
    
    # Remove from cases 
    data.cases <- data.cases %>% .[,c(allele1Rem, allele2Rem)] %>% filter(get(allele1Rem) != alleleRem & get(allele2Rem) != alleleRem)
  }
  
  # Count number of the allele of Interest in data controls
  alleleCount <- alleleCount + table(data.controls[,c(allele1)])[allele] %>% unname() +
    table(data.controls[,c(allele2)])[allele] %>% unname()

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
  alleles <- list(data[, A1], data[, A2]) %>% unlist() 
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  ACFREQ.df <- data.frame(allele = alleles.count %>% names(), 
                           alleleCount = alleles.count %>% as.vector(),
                           alleleFreq = alleles.freq %>% as.vector(), 
                           alleleTotal = nrow(data))
  
  # Remove allele of interest
  ACFREQ.df <- ACFREQ.df %>% filter(allele != settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[2])
  
  # Remove alleles to remove 
  for (A in settings$allele2Remove){
    # Locus and allele
    locusREM <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    c(allele1Rem, allele2Rem) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
    alleleRem = settings$allele2Remove %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Remove 
    ACFREQ.df <- ACFREQ.df %>% filter(allele != alleleRem)
    
  }
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',3), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',3), sep = ''))
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
  fishers.LI.alleles <- c(); fishers.LI.carriers <- c(); OR <- c();
  
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
    #OR
    OR <- c(OR, (cont.table.allele[1]*cont.table.allele[4])/(cont.table.allele[2]*cont.table.allele[3]))
    
  }
  stats.data <- data.frame(allele = ACFREQ.df$allele, 
                           FishersAllelePVAL = fishers.pval.alleles, FishersAlleleOR = fishers.OR.alleles, 
                           FishersAlleleLI = fishers.LI.alleles, FishersAlleleUI = fishers.UI.alleles,
                           ChiAllelePVAL = chi2.pval.alleles, ChiAlleleVAL = chi2.OR.alleles, ORAllele = OR)
  stats.data[is.na(stats.data)] <- 0
  
  return(stats.data)
  
}

######### CHI2 HOMOZYGOSITY / HETEROZYGOSITY HLA ANALYSIS #######

diseases <- c('LGI1')
locus <- paste('HLA', settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[1], sep = '')

for (DS in diseases){
  
  # Get patient information
  system(paste('Rscript --vanilla getPatInfo.R', DS, sep = ' '))
  
  # Pre-preprocess data 
  c(data.cases, data.controls, data.filt, alleleCount) %<-% preprocessData(data, settings) %>% unlist() %>% invisible()
  
  # Compute allele, carrier, and homo/hetero count in cases and controls 
  ACFREQ.cases <- computeACFREQ(data = data.cases, locus = locus, Dx = 'case' ); totalCases <-unique(ACFREQ.cases$alleleTotalCase)
  ACFREQ.controls <- computeACFREQ(data = data.controls, locus = locus, Dx = 'control'); totalControls <- unique(ACFREQ.controls$alleleTotalControl)*2
  ACFREQ.df <- merge(ACFREQ.cases, ACFREQ.controls, by = "allele", all = TRUE); ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df));ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df));
  
  # Compute Chi2
  stats.data <- computeChi2(ACFREQ.df = ACFREQ.df)
  
  # Create dataframes 
  HLA.alleles.df <- merge(stats.data[,c(1,which(grepl('Allele', colnames(stats.data))))],
                          ACFREQ.df[,c(which(grepl(paste(c('allele'), collapse = '|'), colnames(ACFREQ.df))))],
                          by = 'allele')
  
  # Write to excel output
  if (length(settings$allele2Remove) == 0) {
    write.xlsx(x = HLA.alleles.df, file = paste(settings$directory$HLA_Output_HomoHetero,
                                                '/HLA_HardyWeinberg_', DS, '_', settings$alleleOfInterest,'.xlsx', sep = ''),
               col.names = TRUE, row.names = FALSE, append = TRUE)
  } else{
    write.xlsx(x = HLA.alleles.df, file = paste(settings$directory$HLA_Output_HomoHetero,
                                                '/HLA_HardyWeinberg_', DS, '_', settings$alleleOfInterest, '_ControlFor_', settings$allele2Remove, 
                                                '.xlsx', sep = ''),
               col.names = TRUE, row.names = FALSE, append = TRUE)
  }

}
