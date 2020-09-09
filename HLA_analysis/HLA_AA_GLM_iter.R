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
library(plyr)

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

controlAA = function(settings, data, AA_alignment, AA2control){
  
  # Initialize loop 
  GWASID.total <- c(); data.AAControl <- data.frame(GWASID = data$GWASID); 
  if(length(AA2control)>0){
    for (idx in 1:length(AA2control)){
      
      # Get locus, position and AA to control 
      AA <- AA2control[idx] %>% strsplit('_') %>% unlist()
      L2control <- AA[1]; posControl <- AA[2]; AAcontrol <- AA[3]
      c(L2controlA1, L2controlA2) %<-% c(paste0('HLA',L2control,'_A1'), paste0('HLA',L2control,'_A2'))
      
      # Get sequences
      AA2control_locus <- AA_alignment %>% filter(locus == AA[1])
      
      # Get alleles 
      allelesAA.OG <- AA2control_locus %>% apply(MARGIN = 1, function(x, posControl, AAcontrol)
        if (substr(x[3],posControl,posControl) == AAcontrol) {return(x[2])}, posControl, AAcontrol) %>% unlist()
      allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                             .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist()
      
      # Get presence or absence of the allele, append to matrix
      evalPres = function(x){
        if (x[[L2controlA1]] %in% allelesAA | x[[L2controlA2]] %in% allelesAA){
          return(1)
        }else{
          return(0)
        }
      }
      data.AApres <- apply(X = data, MARGIN = 1, FUN = evalPres)
      data.AAControl <- cbind(data.AAControl, data.AApres)
    }
    
    # Rename dataframe columns
    colnames(data.AAControl)[2:ncol(data.AAControl)] <- AA2control
    return(data.AAControl)
  } 
  else{
    return(data.AAControl <- data.frame(GWASID = data$GWASID))
  }
 
}

######### DATA TO ONE HOT ENCODING #########

data2OHE = function(settings, data, allelesAA, PCs, AA_alignment, AA2control){
  
  # Get data with AAs
  data.AApos <- data %>% filter(get(A1) %in% allelesAA| get(A2) %in% allelesAA)
  data.AAneg <- data %>% filter(GWASID %notin% data.AApos$GWASID)
  
  # Create dataframe with OHE
  data.AA <- data.frame(GWASID = c(data.AApos$GWASID, data.AAneg$GWASID), AA = c(rep(1, nrow(data.AApos)), rep(0, nrow(data.AAneg))),
                        Dx = c(data.AApos$Dx, data.AAneg$Dx))
  
  # Merge with Principal components
  data.AA <- merge(data.AA, PCs, by.x = "GWASID", by.y = "ID")
  
  # Merge with presence of controlled amino acid
  data.AAcontrol <- controlAA(settings, data, AA_alignment, AA2control)
  data.AA <- merge(data.AA, data.AAcontrol, by = "GWASID")
  
  # Return 
  return(data.AA)
}

######### RUN GENERALIZED LINEAR MODEL #########

runGLM = function(settings, data.AA, AA2control){
  
  # GLM formula 
  AA.name <- colnames(data.AA)[2]
  if (length(AA2control)==0){
    glm.formula <- paste0("Dx ~ ", paste(AA.name, "PC1","PC2","PC3", sep = ' + '))
  } else{
    glm.formula <- paste0("Dx ~ ", paste(c(AA.name, "PC1","PC2","PC3",AA2control), collapse = ' + '))
    
  }
  
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

# Initialize while lopp
pvalMin <- 0; 
AA.model.df.SIGN <- data.frame(); AA2control <- c()

# While there is a significant Amino acid
idx <- 1
while (pvalMin <= settings$pvalThres){
  print(paste('Current loop:', as.character(idx)), sep = '')
  
  # Initialize dataframe
  AA.model.df <- data.frame()
  
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
          if (paste(L,pos,AA, sep='_') %in% AA2control){
            next
          }
          
          # Get alleles with AAs
          allelesAA.OG <- AA_locus %>% apply(MARGIN = 1, function(x, pos, AA) if (substr(x[3],pos,pos) == AA) {return(x[2])}, pos, AA) %>% unlist()
          allelesAA <- allelesAA.OG %>% lapply(function(x) strsplit(x, split='\\*') %>% unlist() %>% 
                                                 .[2] %>% strsplit(split=':') %>% unlist() %>% .[1:2] %>% paste(collapse=':')) %>% unlist() 
          
          # Crete one hot encoding dataframes
          data.AA <- data2OHE(settings, data, allelesAA, PCs, AA_alignment, AA2control); colnames(data.AA)[2] <- paste(L,pos,AA, sep = '_')
          
          # Run generalized linear model 
          model.AA <- runGLM(settings, data.AA, AA2control)
          
          # If rank is not equal to model matrix (dummy variable trap) skip
          if (is.null(model.AA)){
            next 
          }
          
          # Else, get coefficients
          coefs <- model.AA$coefficients; pvalDim <- dim(coefs)[2]; nVars <- dim(coefs)[1] 
          
          # Append to dataframe 
          AA.model.df<- rbind(AA.model.df, c(L, pos, AA, coefs[2,1], coefs[1:nVars,pvalDim], paste(allelesAA.OG, collapse = ', ')))
        }
      }
    }
  }
  
  # Change colnames
  AA.model.df <- cbind(AA.model.df[,1:6], data.frame(AA.PVALBONF = as.numeric(AA.model.df[,6])*nrow(AA.model.df))) %>%
    cbind(AA.model.df[7:ncol(AA.model.df)])
  colnames(AA.model.df) <- c('Locus', 'POS', 'AA', 'AA.estimate', 'Intercept.PVAL', 'AA.PVAL',"AA.PVALBONF",
                             paste0(rownames(coefs)[3:nVars], rep('.PVAL', nVars-3)),"alleles")
  
  # Check minimum pvalue and append next AA to control
  pvalMin <- min(AA.model.df$AA.PVAL %>% as.numeric()); pvalIdx <- which(AA.model.df$AA.PVAL %>% as.numeric() == pvalMin)[1]
  AA2control <- c(AA2control, paste0(AA.model.df$Locus[pvalIdx],'_', AA.model.df$POS[pvalIdx], '_', AA.model.df$AA[pvalIdx]))
  
  # Append to dataframe of significant alleles
  if (nrow(AA.model.df.SIGN)== 0){
    AA.model.df.SIGN <- rbind(AA.model.df.SIGN, AA.model.df[pvalIdx,])
  } else{
    AA.model.df.SIGN <- rbind.fill(AA.model.df.SIGN, AA.model.df[pvalIdx,])
  }
  
  # Update idx
  idx <- idx +1
  
}

# Save dataframe 
for (L in unique(AA.model.df$Locus)){
  out <- AA.model.df %>% filter(Locus == L)
  write.xlsx(x = out, file = paste(settings$directory$HLA_Output_AA_Analysis, 'HLA_AA_GLM_iter.xlsx', sep = ''), sheetName = L, 
             col.names = TRUE, row.names = FALSE, append = TRUE)
}
write.xlsx(x = AA.model.df.SIGN, file = paste(settings$directory$HLA_Output_AA_Analysis, 'HLA_AA_GLM_iter.xlsx', sep = ''), sheetName = 'significant_AA' , 
           col.names = TRUE, row.names = FALSE, append = TRUE)
write.csv(data.AA, paste(settings$directory$HLA_Output_AA_Analysis, 'dataAA.csv', sep = ''), row.names = FALSE)


########## LOGISTIC GLM PREDICTION ##########

# Formula
glm.formula <- "Dx ~ DRB1_33_Q + DRB1_76_F + DPB1_105_V + DPB1_85_A + DQB1_98_E + B_93_A + 
  DRB4_164_S + C_121_W + DQB1_89_S + A_133_F + DQB1_23_A + DRB1_59_Y + 
  DQB1_58_G + A_119_L + DRB1_210_T + DPB1_62_Q + C_243_R + C_180_R"

glm.formula <- "Dx ~ DRB1_33_Q + DPB1_105_V + DRB1_59_Y  -1"

avgACC <- c(); keepModel <- NULL; maxAcc <- 0
# Data preparation 
for (i in 1:1000){
  data.RAND <- data.AA[sample(1:nrow(data.AA)),]
  data.CASE <- data.RAND %>% filter(Dx == 1); data.CTRL <- data.RAND %>% filter(Dx == 0) %>% .[1:nrow(data.CASE),]
  data.TRAIN <- rbind(data.CASE[1:round(0.9*nrow(data.CASE)),], data.CTRL[1:round(0.9*nrow(data.CTRL)),])
  data.TEST <- rbind(data.CASE[(round(0.9*nrow(data.CASE))+1):nrow(data.CASE),], data.CTRL[(round(0.9*nrow(data.CTRL))+1):nrow(data.CTRL),])
  
  # Train model 
  model.GLM <- glm(data = data.TRAIN, formula =  glm.formula, family = 'binomial', maxit = 100)
  
  # Predict 
  Dx.pred <- predict(model.GLM, newdata = data.TEST, type = 'response')
  
  # Add to test
  data.TESTRES <-  data.frame(GWASID = data.TEST$GWASID, model_prob = Dx.pred, model_pred = 1*(Dx.pred %>% round()), 
                              Dx_test = data.TEST$Dx, good_pred = Dx.pred %>% round() == data.TEST$Dx)
  
  # Calculate accuracy 
  Acc <- sum((1*data.TESTRES$good_pred))/nrow(data.TESTRES) * 100
  
  # Keep model
  if (Acc > maxAcc){
    keepModel <- model.GLM
    keepModelSUM <- keepModel %>% summary()
  }
  
  # Append 
  avgACC <- c(avgACC, Acc)
  
  # Max ACC
  maxAcc <- max(avgACC)
  
}

meanAVG <- mean(avgACC)
keepModelSUM
meanAVG

# Test best model 
Dx.predFINAL <- - predict(keepModel, newdata = data.TEST, type = 'response')

for (i in 1:1000){
  data.RAND <- data.AA[sample(1:nrow(data.AA)),]
  data.CASE <- data.RAND %>% filter(Dx == 1); data.CTRL <- data.RAND %>% filter(Dx == 0) %>% .[1:nrow(data.CASE),]
  data.TRAIN <- rbind(data.CASE[1:round(0.9*nrow(data.CASE)),], data.CTRL[1:round(0.9*nrow(data.CTRL)),])
  data.TEST <- rbind(data.CASE[(round(0.9*nrow(data.CASE))+1):nrow(data.CASE),], data.CTRL[(round(0.9*nrow(data.CTRL))+1):nrow(data.CTRL),])
  
}

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




