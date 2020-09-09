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
library(plyr)
library(TreeBH)

## DESCRIPTION ##
# This scripts uses one-hot encoding to predict with a logistic regression the 
# most significant alleles in each locus in patients where sample size is > 45 

########## IMPORT #########

# Set working directory 
setwd('~/Documents/paraneoplasticSyndromes')

# Import settings 
settings <- jsonlite::fromJSON('HLA_analysis/settingsHLA.json')

# Set values to 0
settings$controlAlleles <- c()
settings$excludePosSubjsByAllele <- c()
settings$allele2Remove <- c()

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
  
  # Carrier frequency
  A1 <- locus %>% paste('_A1', sep = ''); A2 <- locus %>% paste('_A2', sep = '')
  alleles <- list(data[, A1], 
                  data[, A2]) %>% unlist()
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
  ACFREQ.df <- merge(HH.data, carrier.df, by = 'allele')
  
  switch (Dx,
          'case' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Case',4), sep = ''))
            return(ACFREQ.df)
          },
          'control' = {
            colnames(ACFREQ.df) <- c('allele', paste(c('A0','A1','A2','carrierCount', 'carrierFreq', 'carrierTotal'),
                                                     rep('Control',4), sep = ''))
            return(ACFREQ.df)
          }
  )
}

########### CONTROL ALLELES ##########

controlAllele = function(As2control, data.filt){
  
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


runLogisticRegression = function(locus, OHE.carrierFreq.data, PCs, as2control = NULL){
  
  ## Allele Frequency 
  # Merge dataset to include PCs
  PCs$ID <- as.factor(PCs$ID)
  alleles.freq <- colnames(OHE.carrierFreq.data)[-c(1,(ncol(OHE.carrierFreq.data)-length(as2control)):ncol(OHE.carrierFreq.data))]
  OHE.carrierFreq.data <- merge(OHE.carrierFreq.data, PCs, by = 'ID')
  
  # Remove alleles for control 
  locus.subset <- as2control[grepl(as2control, pattern = strsplit(locus, 'HLA') %>% unlist() %>% .[2])]
  allele.subset <- locus.subset %>% lapply(function(x) x %>% strsplit('\\*') %>% unlist() %>% .[2]) %>% unlist()
  OHE.carrierFreq.data[allele.subset] <- NULL
  alleles.freq <- alleles.freq[!alleles.freq %in% allele.subset]
  
  # Run logistic regression on carrier frequency 
  Acarrier.model.df <- data.frame()
  for (allele in alleles.freq){
    if (!is.null(as2control)){
      control.alleles <- paste(' ', as2control %>% sapply(function (x) paste('`', x ,'`', sep = '')) %>% paste(collapse = ' + '), sep = '+ ')
    } else{
      control.alleles <- ''
    }
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
  glm.data <- Acarrier.model.df
  
  # Return
  return(glm.data)
  
}

fitGLM = function(settings, locus, data.filt, data.cases, data.controls, PCs, as2control = NULL){
  
  # Control for allele
  controlAllele.df <- controlAllele(as2control, data.filt)
  
  # Subset locus
  allele1 <- paste(locus, '_A1',sep = '')
  allele2 <- paste(locus,'_A2', sep =  '')
  data.locus <- data.filt[,c('GWASID', allele1, allele2)]
  
  # Compute allele frequencies and coutns, and carrier frequencies and counts, and merge
  # Compute allele and carrier count and frequency, and merge
  ACFREQ.cases <- computeACFREQ(data.cases, locus, 'case');
  carrierCases <- unique(ACFREQ.cases$carrierTotalCase)
  ACFREQ.controls <- computeACFREQ(data.controls, locus, 'control');
  carrierControls <- unique(ACFREQ.controls$carrierTotalControl)
  ACFREQ.df <- merge(ACFREQ.cases[,!names(ACFREQ.cases) %in% c('A0','A1','A2')], ACFREQ.controls, by = 'allele', all = TRUE) 
  ACFREQ.df[is.na(ACFREQ.df)] <- 0
  ACFREQ.df$carrierTotalCase <- rep(carrierCases, nrow(ACFREQ.df)); ACFREQ.df$carrierTotalControl <- rep(carrierControls, nrow(ACFREQ.df))
  ACFREQ.df <- ACFREQ.df %>% filter(allele != '')
  
  # Parse one hot encoding and merge
  OHE.carrierFreq.data <- carrierFreqOHE(data.locus); if ('V1' %in% colnames(OHE.carrierFreq.data)){OHE.carrierFreq.data$V1 <- NULL}
  OHE.carrierFreq.data<-  merge(as.data.frame(OHE.carrierFreq.data), data.filt[c('Dx', 'GWASID')], by.x ='ID', by.y = 'GWASID') %>% 
    merge(controlAllele.df, by = 'ID')
  
  # Remove subjects thar are not controls or cases
  OHE.carrierFreq.data <- OHE.carrierFreq.data %>% filter(Dx != -9)
  
  # Run logistic regression model 
  glm.data <- runLogisticRegression(locus, OHE.carrierFreq.data, PCs, as2control)
  
  # Create dataframes 
  HLA.GLM_carriers.df <-merge(glm.data[,c(1,which(grepl('CARRIER', colnames(glm.data))))],
                              ACFREQ.df[,c(1,which(grepl(paste(c('A0','A1','A2', 'carrier'), collapse = '|'), colnames(ACFREQ.df))))],
                              by = 'allele')
  
  return(HLA.GLM_carriers.df)
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
  
  # Initialize while lopp
  pval <- 0; as2control <- c(); signAlleles <- list(); 
  
  #### First pass 
  # Get locus of interest
  locusOfInterest <- settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[1]
  locusOfInterest <- sapply(locusOfInterest, function(x) paste('HLA',x, sep=''))
  
  # For each locuts fit a GLM
  for (locus in locusOfInterest){
        HLA.GLM_carriers.df <- fitGLM(settings, locus, data.filt, data.cases, data.controls, PCs)
  }
  
  # Add to outputs
  as2control <- c(settings$alleleOfInterest)
  preOut <- HLA.GLM_carriers.df %>% filter(allele == settings$alleleOfInterest %>% strsplit('\\*') %>% unlist() %>% .[2])
  preOut$allele <- settings$alleleOfInterest; signAlleles[[1]] <- preOut
  
  # HLA Loci
  loci <- paste(c(rep('HLA', 10)), c('A','B','C','DQA1', 'DQB1', 'DPB1', 'DRB1','DRB3','DRB4','DRB5'), sep = '')
  
  # While signifiant alleles
  idx <- 2; iter <- 1
  while (pval < 0.05){
    
    # Fit GLM 
    HLA.GLM_carriers.list <- list(); pvalTotal <- c()
    for (locus in loci){
      HLA.GLM_carriers.df <- fitGLM(settings, locus, data.filt, data.cases, data.controls, PCs, as2control)
      HLA.GLM_carriers.list[[locus]] <- HLA.GLM_carriers.df 
      for (A in HLA.GLM_carriers.df$allele){
        pvalA <- HLA.GLM_carriers.df$allele.CARRIER.pval[which(HLA.GLM_carriers.df$allele == A)]
        pvalTotal[paste(str_split(locus, 'HLA') %>% unlist() %>% .[2], '*', A, sep = '')] <- pvalA
      }
    }
    
    # Get minimum allele value 
    pval <- pvalTotal[which(pvalTotal == min(pvalTotal))][1]
    pvalMin <- pval %>% names(); pvalMinLocus <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[1]
    pvalMinAllele <- pvalMin %>% strsplit('\\*') %>% unlist() %>% .[2]
    
    # Add to outputs
    as2control <- c(as2control, pvalMin)
    preOut <- HLA.GLM_carriers.list[[paste('HLA', pvalMinLocus, sep = '')]] %>% filter(allele == pvalMinAllele)
    preOut$allele <- pvalMin  ; signAlleles[[idx]] <- preOut
    
    # Compute allele and locus group for correction
    alleleGroup <- pvalTotal %>% names(); locusGroup <- lapply(alleleGroup, function(x) x %>% strsplit('\\*') %>% unlist() %>% .[1]) %>% unlist()
    groups <- matrix(c(as.factor(locusGroup), as.factor(alleleGroup)), ncol = 2)
    
    # Compute TreeBH
    TreeBH.res <- get_TreeBH_selections(pvals = pvalTotal %>% unname(), groups = groups, q = rep(0.05, ncol(groups)))
    TreeBH.res.DF <- data.frame(allele = alleleGroup, pval = pvalTotal %>% unname()); 
    TreeBH.res.DF <- TreeBH.res.DF[which(TreeBH.res[,ncol(TreeBH.res)] == 1),]
    
    # Write TreeBH for each iteration, and append at the end of the sheet controlled alleles
    if (nrow(TreeBH.res.DF) > 0){
      write.xlsx(x = TreeBH.res.DF,  paste(settings$directory$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter_', DS,'.xlsx', sep = ''), 
                 sheetName = paste('Iteration', as.character(iter), sep = ''), col.names = TRUE, row.names = FALSE )
      write.xlsx(x = '\n', paste(settings$directory$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter_', DS,'.xlsx', sep = ''), 
                 sheetName = paste('Iteration', as.character(iter), sep = ''), append = TRUE);
      for (A in as2control){
        write.xlsx(x = as.character(A), paste(settings$directory$HLA_Output_GLM_Iter, 'HLA_TreeBH_iter_', DS,'.xlsx', sep = ''), 
                   sheetName = paste('Iteration', as.character(iter), sep = ''), append = TRUE);
      }
    }
    
    # Update
    iter <- iter + 1
    idx <- idx +1
    
  }
    
    # Format output
    allele <- as.data.frame(signAlleles[[1]]); 
    for (idx in 2:length(signAlleles)){
      allele <- rbind.fill(allele, signAlleles[[idx]] %>% as.data.frame())
    }
    
    # Write
    for (idx in 1:length(HLA.GLM_carriers.list)){
      HLA.GLM_carriers.df <- HLA.GLM_carriers.list[[idx]]; locus <- HLA.GLM_carriers.list %>% names() %>% .[idx]
      write.xlsx(x = HLA.GLM_carriers.df, file = paste(settings$directory$HLA_Output_GLM_Iter, 'HLA_GLM_Carriers_', DS,'.xlsx', sep = ''), sheetName = locus, 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
    }
    write.xlsx(x = allele, file = paste(settings$directory$HLA_Output_GLM_Iter, 'HLA_GLM_Carriers_', DS,'.xlsx', sep = ''), sheetName = 'Significant_alleles', 
               col.names = TRUE, row.names = FALSE, append = TRUE)
  
}

