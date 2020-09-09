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
# This scripts uses one-hot encoding to predict with a logistic regression the 
# most significant alleles in each locus in patients where sample size is > 45 

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


# Import Principal Components 
if (settings$allPlates == 1){
  # Import provided PCs 
  PCs <- read.table(settings$file$PCS_all_plates, header = TRUE)
  PCs$ID <- PCs$IID; PCs <- select(PCs, -c('IID', 'FID'))
  changeIdFormat = function(GWASID){
    if (grepl('Stanford', GWASID)){
      idpre <- GWASID %>% strsplit('_')  %>% unlist() 
      return(strsplit(idpre[length(idpre)], '\\.')[[1]][1])
    } else{
      return (GWASID)
    }
  }
  PCs$ID <- sapply(PCs$ID, changeIdFormat) %>% unlist()
  
} else{
  PCs <- read.table(file = settings$file$PCS, header = TRUE)
  PCs <- cbind(ID = sapply(PCs$IID %>% as.character(), function(x) strsplit(x, '_') %>% unlist() 
                     %>% tail(1) %>% strsplit('\\.') %>% unlist() %>% head(1)),
               PCs, row.names= NULL)
  PCs <- select(PCs, -c('IID', 'FID'))
}


########### ONE HOT ENCODING FUNCTIONS ###############

# HLA Parse function
alleleFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  ID = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(ID, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, ID~make_HLA, fun.aggregate = length)
  dcast_HLA$ID = as.factor(dcast_HLA$ID)
  return(setDF(dcast_HLA))
}

# HLA Parse function
carrierFreqOHE=function(test_DF){
  test_DF = as.data.frame(test_DF)
  ID = rep(test_DF[,1], 2)
  make_HLA = c(as.character(test_DF[,2]), as.character(test_DF[,3]))
  make_DF = cbind.data.frame(ID, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, ID~make_HLA, fun.aggregate = length)
  dcast_vals <- dcast_HLA[,!c("ID")]; dcast_vals[dcast_vals > 1] <- 1
  dcast_HLA <- cbind(data.frame(ID =dcast_HLA$ID), dcast_vals)
  dcast_HLA$ID = as.factor(dcast_HLA$ID)
  return(setDF(dcast_HLA))
}

########### PRE-PROCES DATA ############
# Compute a logistic regression to test significance of carrier frequencies 

## Exclude subjects w/ positive allele 
excludeByAllele = function(data.filt, settings){
  
  # Alleles to control
  As2exclude = settings$excludePosSubjsByAllele
  
  exclude_list <- c()
  for (A in As2exclude){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    A1 <- paste('HLA', locus, '_A1', sep = ''); A2 <- paste('HLA', locus, '_A2', sep = '')
    allele2exclude = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Get data 
    data.filt <- data.filt %>% filter(get(A1) != allele2exclude & get(A2) != allele2exclude)
    
  }
  return(data.filt)
}

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
  
  # Exclude subjects with allele
  data.filt <- excludeByAllele(data.filt, settings)
  
  # Return
  return(data.filt)
  
}

########### ALLELE AND CARRIER FREQUENCY  #################


computeACFREQ = function(data, locus, Dx){
  # Compute heterozigous, homozigous, or absence count. 
  # Compute allele frequency, count and total.
  # Compute carrier frequency, count and total.
  
  # Allele frequency 
  A1 <- locus %>% paste('_A1', sep = ''); A2 <- locus %>% paste('_A2', sep = '')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
  alleles.count <- table(alleles)
  alleles.freq <- alleles.count %>% prop.table() * 100
  alleles.df <- data.frame(allele = alleles.count %>% names(), 
                           alleleCount = alleles.count %>% as.vector(),
                           alleleFreq = alleles.freq %>% as.vector(), 
                           alleleTotal = nrow(data)*2)
  
  # Carrier frequency
  carriers <- data[,c(A1, A2)]
  carriers.levels <- list(data[, A1], 
                          data[, A2]) %>% unlist() %>% levels()
  carriers.unique <- apply(carriers, 1, function(x) unique(x)) %>% unlist() %>% as.factor()
  carriers.count <- table (carriers.unique); carriers.count[c(carriers.levels %>% setdiff(carriers.count %>% names()))] <- 0
  carriers.freq <- carriers.count /nrow(data) * 100
  carrier.df <- data.frame(allele = carriers.count %>% names(), 
                           carrierCount = carriers.count %>% as.vector(),
                           carrierFreq = carriers.freq %>% as.vector(),
                           carrierTotal = nrow(data))
  
  # Heterozigous, homozigous, and absence count
  A0 <- c(); A1 <- c(); A2 <- c();
  for (A in levels(as.factor(alleles))){
    HH.data <- carriers[which(carriers[,1]==as.character(A) | carriers[,2]==as.character(A)),]
    HH.count <- HH.data %>% apply(1, function(x) x %>% unlist() %>% unique() %>% length()) %>% unlist() %>% table()
    A0 <- c(A0, nrow(data) - nrow(HH.data)); A1 <- c(A1,HH.count['2'] %>% unname()); A2 <- c(A2,HH.count['1'] %>% unname())
  }
  A1[is.na(A1)] <- 0; A2[is.na(A2)] <- 0; HH.data <- data.frame(allele = levels(as.factor(alleles)), A0 = A0, A1 = A1, A2 = A2)
  
  # Merge 
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele') %>% merge(alleles.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Case',7), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal',
                                                       'alleleCount', 'alleleFreq','alleleTotal'),
                                                     rep('Control',7), sep = ''))
            return(ACFREQ.df)
          }
  )
}

########### CONTROL ALLELES ##########

controlAllele = function(settings, locus, data.filt){
  
  # Alleles to control
  As2control = settings$controlAlleles
  
  # Get unique loci 
  lociControl = lapply(As2control, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist() %>% unique() 
  
  # For each allele 
  alleleControl.df = data.frame(ID = data.filt$GWASID)
  for (A in As2control){
    
    # Get locus and allele 
    locus <- A %>% strsplit('\\*') %>% unlist() %>% .[1]
    allele2control = A %>% strsplit('\\*') %>% unlist() %>% .[2]
    allele1 <- paste('HLA', locus, '_A1', sep = '')
    allele2 <- paste('HLA', locus, '_A2', sep = '')
    
    # Get subjects
    alleleControl.df[A] <- as.logical(c(data.filt[,c(allele1)] %>% as.character() == allele2control) + 
                            c(data.filt[,c(allele2)] %>% as.character() == allele2control)) %>% as.integer()
    
  }
  return(alleleControl.df)
}


############ REGRESSION MODEL ############


runLogisticRegression = function(locus, OHE.alleleFreq.data, OHE.carrierFreq.data, PCs){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  as2control <- settings$controlAlleles
  alleles.freq <- colnames(OHE.alleleFreq.data)[-c(1,(ncol(OHE.alleleFreq.data)-length(as2control)):ncol(OHE.alleleFreq.data))]
  OHE.alleleFreq.data <- merge(OHE.alleleFreq.data, PCs, by = 'ID')
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, PCs, by = 'ID')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = strsplit(locus, 'HLA') %>% unlist() %>% .[2])]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.alleleFreq.data[allele.subset] <- NULL
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on allele frequency data
  Afreq.model.df <- data.frame()
  for (allele in alleles.freq){
    control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    glm.formula <- paste('Dx ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    Afreq.model <- glm(data = OHE.alleleFreq.data, 
                       formula = as.formula(glm.formula),
                       family = 'binomial', maxit = 100) %>% summary()
    Afreq.model.df <- rbind(Afreq.model.df, c(Afreq.model$coefficients[2,1], 
                                              Afreq.model$coefficients[,dim(Afreq.model$coefficients)[2]]))
  }
  colnames(Afreq.model.df) <- c('allele.COEF.ALLELE', 
                                c('Incercept', 'allele', Afreq.model$coefficients[-c(1,2),] %>% row.names()) %>% paste('.ALLELE.pval', sep = ''))
  Afreq.model.df <- cbind(data.frame(allele=alleles.freq, Afreq.model.df))
  
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  for (allele in alleles.freq){
    control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    glm.formula <- paste('Dx ~ `',allele, '` + PC1 + PC2 + PC3', control.alleles, sep = '')
    Acarrier.model <- glm(data = OHE.carrierFreq.data, 
                       formula = as.formula(glm.formula),
                       family = 'binomial', maxit = 100) %>% summary()
    Acarrier.model.df <- rbind(Acarrier.model.df, c(Acarrier.model$coefficients[2,1], 
                                                    Acarrier.model$coefficients[,dim(Acarrier.model$coefficients)[2]]))
  }
  colnames(Acarrier.model.df) <- c('allele.COEF.CARRIER', c('Incercept', 'allele', Acarrier.model$coefficients[-c(1,2),] %>% row.names()) %>%
                                  paste('.CARRIER.pval', sep = ''))
  Acarrier.model.df <- data.frame(allele=alleles.freq, Acarrier.model.df)
  
  # Merge 
  glm.data <- merge(Acarrier.model.df, Afreq.model.df, by = 'allele')
  
  # Return
  return(glm.data)
  
}

############### HLA ANALYSIS ###########

# Diseases 
diseases <- c('LGI1')

for (DS in diseases){
  
  # Get patient data
  system(paste('Rscript --vanilla getPatInfo.R', DS, sep = ' '))
  
  # Import controls and cases, and merge
  if (settings$allPlates != 1) {
    cases <- read.table(settings$file$GWASIDsCases) %>% unlist()
    controls <- read.table(settings$file$GWASIDsControls) %>% unlist()
    
    # Pre-process data
    data.filt <- preprocessData(cases, controls)
  } else {
    
    # Data already pre-processed
    data.filt <- excludeByAllele(data, settings)
  }
  
  # Get CASES and CONTROLS 
  data.cases <- data.filt %>% filter(Dx == 1)
  data.controls <- data.filt %>% filter(Dx == 0)
  
  # HLA Loci
  loci <- paste(c(rep('HLA', 10)), c('A','B','C','DQA1', 'DQB1', 'DPB1', 'DRB1','DRB3','DRB4','DRB5'), sep = '')
  
  # Control for allele
  controlAllele.df <- controlAllele(settings, loci, data.filt)
  
  # Iterate over loci for univariate analysis
  models.df <- data.frame()
  for (locus in loci){
    
    # Subset locus
    allele1 <- paste(locus, '_A1',sep = '')
    allele2 <- paste(locus,'_A2', sep =  '')
    data.locus <- data.filt[,c('GWASID', allele1, allele2)]
    
    # Compute allele frequencies and coutns, and carrier frequencies and counts, and merge
    # Compute allele and carrier count and frequency, and merge
    ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
    totalCases <-unique(ACFREQ.cases$alleleTotalCase); carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
    ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
    totalControls <-unique(ACFREQ.controls$alleleTotalControl); carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
    ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
    ACFREQ.df[is.na(ACFREQ.df)] <- 0
    ACFREQ.df$alleleTotalCase <- rep(totalCases, nrow(ACFREQ.df)); ACFREQ.df$alleleTotalControl <- rep(totalControls, nrow(ACFREQ.df));
    ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
    ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
    
    # Parse one hot encoding and merge
    OHE.alleleFreq.data <- alleleFreqOHE(data.locus)
    OHE.alleleFreq.data <- merge(as.data.frame(OHE.alleleFreq.data), data.filt[c('Dx', 'GWASID')], by.x ='ID', by.y = 'GWASID') %>% 
      merge(controlAllele.df, by = 'ID')
    OHE.carrierFreq.data <- carrierFreqOHE(data.locus)
    OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), data.filt[c('Dx', 'GWASID')], by.x ='ID', by.y = 'GWASID') %>% 
      merge(controlAllele.df, by = 'ID')
    
    # Remove subjects thar are not controls or cases
    OHE.alleleFreq.data <- OHE.alleleFreq.data %>% filter(Dx != -9)
    OHE.carrierFreq.data <- OHE.carrierFreq.data %>% filter(Dx != -9)
    
    # Run logistic regression model 
    glm.data <- runLogisticRegression(locus, OHE.alleleFreq.data, OHE.carrierFreq.data, PCs)
    
    # Create dataframes 
    HLA.GLM_alleles.df <-merge(glm.data[,c(1,which(grepl('ALLELE', colnames(glm.data))))],
                                ACFREQ.df[,c(which(grepl(paste(c('A0','A1','A2', 'allele'), collapse = '|'), colnames(ACFREQ.df))))],
                                by = 'allele')
    HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                                ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                                by = 'allele')
    
    # Write to excel output
    if (length(settings$excludePosSubjsByAllele) > 0){
      write.xlsx(x = HLA.GLM_alleles.df, file = paste(settings$directory$HLA_Output_GLM, 'HLA_GLM_Alleles_', DS,'_',
                                                      settings$excludePosSubjsByAllele,'_NEGATIVE','.xlsx', sep = ''), sheetName = locus, 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
      write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$directory$HLA_Output_GLM, 'HLA_GLM_Carriers_',DS,'_',
                                                       settings$excludePosSubjsByAllele,'_NEGATIVE','.xlsx', sep = ''), sheetName = locus, 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
    } else {
    write.xlsx(x = HLA.GLM_alleles.df, file = paste(settings$directory$HLA_Output_GLM, 'HLA_GLM_Alleles_', DS,'.xlsx', sep = ''), sheetName = locus, 
               col.names = TRUE, row.names = FALSE, append = TRUE)
    write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$directory$HLA_Output_GLM, 'HLA_GLM_Carriers_', DS,'.xlsx', sep = ''), sheetName = locus, 
               col.names = TRUE, row.names = FALSE, append = TRUE)
    }
  }
}
