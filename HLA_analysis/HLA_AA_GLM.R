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
library(TreeBH)
library(zeallot)
library(epitools)

########## IMPORT #########

# Operator
`%notin%` <- Negate(`%in%`)

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

# Import amino acid alignment 
AA_alignment <- read.table(settings$file$Output_AminoacidAlignment_PostPro, header = TRUE, sep = ',')

######### CONTROL FOR AMINO ACID ##########

controlAA = function(settings, data, AA_alignment){
  
  # Get locus, position and AA to control 
  AA2control <- settings$AA2control %>% strsplit('_') %>% unlist()
  L2control <- AA2control[1]; posControl <- AA2control[2]; AAcontrol <- AA2control[3]
  c(L2controlA1, L2controlA2) %<-% c(paste0('HLA',L2control,'_A1'), paste0('HLA',L2control,'_A2'))
  
  # Get sequences
  AA2control_locus <- AA_alignment %>% filter(locus == AA2control[1])
  
  # Get alleles 
  allelesAA.OG <- AA2control_locus %>% apply(MARGIN = 1, function(x, posControl, AAcontrol)
    if (substr(x[3],posControl,posControl) == AAcontrol) {return(x[2])}, posControl, AAcontrol) %>% unlist()
  allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                         .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist()
  
  # Get data and count cases and controls
  data.AApos <- data %>% filter(get(L2controlA1) %in% allelesAA| get(L2controlA2) %in% allelesAA)
  data.AAneg <- data %>% filter(GWASID %notin% data.AApos$GWASID)
  
  # Create dataframe with presence of Amino Acid
  data.AAControl <- data.frame(GWASID = c(data.AApos$GWASID, data.AAneg$GWASID),
                               pos = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))))
  colnames(data.AAControl) <- c('GWASID', settings$AA2control)
  
  return(data.AAControl)
}

######### DATA TO ONE HOT ENCODING #########

data2OHE = function(settings, data, allelesAA, PCs, AA_alignment){
  
  # Get data with AAs
  data.AApos <- data %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
  data.AAneg <- data %>% filter(GWASID %notin% data.AApos$GWASID)
  
  # Create dataframe with OHE
  data.AA <- data.frame(GWASID = c(data.AApos$GWASID, data.AAneg$GWASID), AA = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))),
                        Dx = c(data.AApos$Dx, data.AAneg$Dx))
                        
  # Merge with Principal components
  data.AA <- merge(data.AA, PCs, by.x = "GWASID", by.y = "ID")
  
  # Merge with presence of controlled amino acid
  data.AAcontrol <- controlAA(settings, data, AA_alignment)
  data.AA <- merge(data.AA, data.AAcontrol, by = "GWASID")
  
  # Return 
  return(data.AA)
}

######### RUN GENERALIZED LINEAR MODEL #########

runGLM = function(settings, data.AA){
  
  # GLM formula 
  AA.name <- colnames(data.AA)[2]
  glm.formula <- paste0("Dx ~ ", AA.name, " + PC1 + PC2 + PC3 + ", settings$AA2control)
  
  # Run GLM 
  AA.model <- glm(data = data.AA, formula = as.formula(glm.formula), family = 'binomial', maxit = 100) 
  modelMat <- model.matrix(AA.model)
  
  # Check if the model matrix has full and return coefficients, otherwise return NULL
  rank <- qr(modelMat)$rank
  if (rank == ncol(modelMat)){
    return(AA.model %>% summary())
  } else {
    return(NULL)
  }
}

########## AMINO ACID ANALYSIS ##########

# Get loci, number of cases and controls
loci <- unique(AA_alignment$locus)
c(Ncases, Ncontrols) %<-% c(data$Dx %>% table %>%.['1'], data$Dx %>% table %>%.['0'])

# Initialize loop
locus <- c(); AA.eval <- c();  pos.eval <- c(); AA.estimate <- c()
intercept.PVAL <- c(); AA.pval <- c(); PC1.pval <-c(); PC2.pval <- c(); PC3.pval <- c(); AA2control.pval <- c();
alleles <- c()

# For each locus 
for (L in loci){
  
  # Get locus ID
  c(A1, A2) %<-% c(paste0('HLA',L,'_A1'), paste0('HLA',L,'_A2'))
  
  # Get alignment subset, and counts subset
  AA_locus <- AA_alignment %>% filter(locus == L)
  
  # Get max sequence length 
  maxLen <- AA_locus$sequence %>% lapply(nchar) %>% unlist() %>% max()
  
  # Get unique AAs 
  for (pos in 1:maxLen){
    
    # Get unique AAs
    AAs <- AA_locus$sequence %>% lapply(function(x, pos) substr(x, pos, pos), pos) %>% unlist() %>% unique()
    AAs <- AAs[which(AAs != '' & AAs!= '*')]
    
    # If more than one, count 
    if (length(AAs >1)){
      
      for (AA in AAs){
        
        # Skip controlled alleles
        if (paste(L,pos,AA, sep='_') == settings$AA2control){
          next
        }
        
        # Get alleles with AAs
        allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
        allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                               .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 

        # Crete one hot encoding dataframes
        data.AA <- data2OHE(settings, data, allelesAA, PCs, AA_alignment); colnames(data.AA)[2] <- paste(L,pos,AA, sep = '_')
        
        # Run generalized linear model 
        model.AA <- runGLM(settings, data.AA)
        
        # If rank is not equal to model matrix (dummy variable trap) skip
        if (is.null(model.AA)){
          next 
        }
        
        # Else, get coefficients
        coefs <- model.AA$coefficients; pvalDim <- dim(coefs)[2]
        
        # Append to vectors
        locus <- c(locus, L); AA.eval <- c(AA.eval, AA); pos.eval <- c(pos.eval, pos);
        AA.estimate <- c(AA.estimate, coefs[2,1]); intercept.PVAL <- c(intercept.PVAL, coefs[1, pvalDim])
        AA.pval <- c(AA.pval, coefs[2, pvalDim]); PC1.pval <- c(PC1.pval, coefs[3, pvalDim]);
        PC2.pval <- c(PC2.pval, coefs[3, pvalDim]); PC3.pval <- c(PC3.pval, coefs[4, pvalDim]); 
        AA2control.pval <- c(AA2control.pval, coefs[4, pvalDim]);
        alleles <- c(alleles, paste(allelesAA.OG, collapse = '; '))
      }
    }
  }
}

# Create dataframe
AA.analysis.results <- data.frame(locus= locus, AA = AA.eval, pos = pos.eval, AA.estimate = AA.estimate, intercept.PVAL = intercept.PVAL, 
                                  AA.pval = AA.pval, AA.pvalBONF= AA.pval*length(AA.pval), PC1.pval = PC1.pval, PC2.pval = PC2.pval, PC3.pval = PC3.pval, AA2control.pval = AA2control.pval,
                                  alleles = alleles)
colnames(AA.analysis.results)[ncol(AA.analysis.results)-1] = settings$AA2control

########## EXTRA FUNCTIONS ############

carrierCount = function(settings, data){
  
  # Get loci, number of cases, number of controls
  loci <- colnames(data)[which(grepl('HLA', colnames(data)))] %>%
    lapply(function(x) x %>% strsplit('_') %>% unlist() %>% .[1] %>% strsplit('HLA') %>% unlist() %>% .[2]) %>% unlist() %>% unique()
  c(Ncases, Ncontrols) %<-% c(data$Dx %>% table() %>% .['1'], data$Dx %>% table() %>% .['0'])
  
  # Initialize loop 
  alleles <- c(); locusVec <- c(); cases.Count <- c(); controls.Count <- c()
  
  # For each locus, count carrier in cases and controls
  for (locus in loci){
    
    # Get alleles 
    c(A1, A2) %<-% c(paste('HLA', locus, '_A1', sep = ''), paste('HLA', locus, '_A2', sep = ''))
    data.locus <- data[, c('GWASID', A1, A2, 'Dx')]
    
    # Get carriers and count 
    carriers.locus <- c(data.locus[,A1], data.locus[,A2]) %>% unique()
    carriers.locus <- carriers.locus[which(carriers.locus != "")]
    for(carrier in carriers.locus){
      data.carrier <- data.locus %>% filter(get(A1) == carrier | get(A2) == carrier)
      c(C.cases, C.controls ) %<-% c(table(data.carrier$Dx)['1'], table(data.carrier$Dx)['0']);
      if(C.cases %>% is.na()){C.cases <- 0}; if (C.controls %>% is.na()) {C.controls <- 0}
      c(F.cases, F.controls) %<-% c(C.cases/Ncases, C.controls/Ncontrols)
      
      # If frequency is smaller than a threshold, next. Else add to vector 
      if (is.na(C.cases) | is.na(C.controls)){
        next 
      } else if (F.cases < settings$min_carrierFreq | F.controls < settings$min_carrierFreq){
        next
      } else{
        locusVec <- c(locusVec, locus)
        A <- paste( locus, '*', carrier, ':01', sep = ''); alleles <- c(alleles, A)
        cases.Count <- c(cases.Count, C.cases); controls.Count <- c(controls.Count, C.controls)
      }
    }
  }
  
  # Create data frame 
  counts.data <- data.frame(locus = locusVec, allele = alleles, Ncases = cases.Count, Ncontrols = controls.Count)
  
  # Return 
  return (counts.data)
  
}




