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

########## IMPORT ##########
setwd("~/Documents/anti-HU")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

########## COLLAPSE #########
path <- settings$directory$HLA_predictions
pred.files <- list.files(settings$directory$HLA_predictions)

# Load first file and get list of samples
init.data <- get(load(paste(path, pred.files[1], sep = '')))
sample.id <- pred.guess$value$sample.id
data <- data.frame('FID' = sample.id, 'IID' = sample.id)

# Collapse
for (file in pred.files){
  temp <- get(load(paste(path, file, sep = '')))
  name <- strsplit(file, '\\.')[[1]][1]
  data[paste(name, 'allele1', sep = '_')] <- as.factor(temp$value$allele1)
  data[paste(name, 'allele2',sep = '_')] <- as.factor(temp$value$allele2)
}

# Write 
write.table(data, file = settings$file$HLA_total, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

############ MERGE #############

# Read eigenvectors
eigenvec <- read.table(settings$file$finalPCAeigenvectors, header = TRUE, sep = "\t", quote = '')

# Merge 
covars <- merge(eigenvec, data, by = c('FID', 'IID'))
covars.names <- colnames(covars)
covars.names <- paste(covars.names, collapse = ', ')
# Write 
write.table(covars, file = settings$file$covariates, quote = FALSE, sep = '\t', row.names = FALSE)
write(covars.names, file = 'Resources/covars_names.txt')

covars.names


########## PREDICT #############

# Import controls and cases, and merge
cases <- read_csv(settings$file$GWASIDsCases, col_names = FALSE)
colnames(cases) <- c('Dx')
controls <- read_csv(settings$file$GWASIDsControls, col_names = FALSE)
colnames(controls) <- c('Dx')
subjects <- rbind(cases, controls)

# Merge to data 
Dx <- c()
ID <- c()
for (subj in data$FID){
  subj.ID <- strsplit(tail(strsplit(as.character(subj),'\\_')[[1]], n = 1), '\\.')[[1]][1]
  if (subj.ID %in% controls$Dx){
    Dx <- c(Dx, 0)
  } else if (subj.ID %in% cases$Dx) {
    Dx <- c(Dx, 1)
  } else{
    Dx <- c(Dx, -9)
  }
  ID <- c(ID, subj.ID)
}

# Add to datframe 
data['Dx'] <- Dx
data['ID'] <- ID

# Filt data 
ctrls.data <- data[which(data$Dx == 0),]
rows <- sample(nrow(ctrls.data))
ctrls.data <- ctrls.data[rows,]
ctrls.data <- ctrls.data[c(1:nrow(cases)),]
cases.data <- data[which(data$Dx == 1),]
bal.data <- rbind(cases.data, ctrls.data)

# Predict unbalanced dataset
data <- data[which(data$Dx != -9),]
data$Dx <- as.factor(data$Dx)
colnames(data)
preds <- ''
for (name in colnames(data)){
  preds <- paste(preds, name, sep = ' + ')
}
preds
model <- glm(Dx ~  HLAA_allele1 + HLAA_allele1*HLAA_allele2 -FID - IID -Dx -ID,
             data=data, family = 'binomial',  control = list(maxit = 2000))
summary(model)

# Predict with balanced dataset
bal.data$Dx <- as.factor(bal.data$Dx)
model.bal <- glm(Dx ~ HLAA_allele1, 
                 data=bal.data, family = 'binomial',  control = list(maxit = 2000))
summary(model.bal)

# Plot 
p = ggplot(bal.data,aes(x=HLAA_allele1,  
                    y=HLAB_allele2,
                    color=Dx))
p + geom_jitter(alpha=1) +  
  scale_color_manual(breaks = c('1','0'),
                     values=c('darkgreen','red'))

########## Functions ##############

PARSE_HLA=function(test_DF){
  test_DF = as.data.frame(test_DF)
  repID = rep(test_DF[,1], 2)
  make_HLA = c(test_DF[,2], test_DF[,3])
  make_DF = cbind.data.frame(repID, make_HLA)
  setDT(make_DF)
  dcast_HLA = dcast(make_DF, repID~make_HLA, fun.aggregate = length)
  dcast_HLA$repID = as.character(dcast_HLA$repID)
  return(dcast_HLA)
  #return(dcast_HLA)#[, -1])
}s