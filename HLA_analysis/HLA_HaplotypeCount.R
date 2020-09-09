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
library(haplo.stats)
library(questionr)
library(epitools)

######### DESCRIPTION #########
# Compute an heterozygote comparison of the allele of interest (e.g. DRBB1*:07:01 = DR7).
# Take DR7 positive heterozygote cases and controls, and compare the allele count between them.


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

# Get cases and controls data
data.cases <- data %>% filter(Dx == 1)
data.controls <- data %>% filter(Dx == 0)

# Create not in operator 
`%notin%` <- Negate(`%in%`)


########### COUNT HAPLOTYPE ###########

# Compute geno files for controls and cases, and labels
geno <- data[,c('HLADRB1_A1', 'HLADRB1_A2', 'HLADRB4_A1','HLADRB4_A2','HLADQB1_A1','HLADQB1_A2','HLADQA1_A1','HLADQA1_A2')]
geno.case <- data.cases[,c('HLADRB1_A1', 'HLADRB1_A2', 'HLADRB4_A1','HLADRB4_A2','HLADQB1_A1','HLADQB1_A2','HLADQA1_A1','HLADQA1_A2')]
geno.controls <- data.controls[,c('HLADRB1_A1', 'HLADRB1_A2', 'HLADRB4_A1','HLADRB4_A2','HLADQB1_A1','HLADQB1_A2','HLADQA1_A1','HLADQA1_A2')]
labels <- c('HLADRB1','HLADRB4','HLADQB1','HLADQA1')

# Compute haplotype couting
haplo.count <- haplo.em(geno = geno, locus.label = labels, miss.val = c(0,NA))
haplo.count.groups <- haplo.group(group = data$Dx, geno, locus.label = labels, miss.val = c(0,NA))

# Sort by maximum frequency in cases
haplo.Ord <- haplo.count.groups$group.df[order(-haplo.count.groups$group.df$`data$Dx=1`),]
haplo.Ord$`data$Dx=0` <- haplo.Ord$`data$Dx=0`*100; haplo.Ord$`data$Dx=1` <- haplo.Ord$`data$Dx=1` *100
haplo.Ord$Total <- haplo.Ord$Total *100
haplo.Ord$`data$Dx=0` <- round(haplo.Ord$`data$Dx=0`, digits = 8); haplo.Ord$Total <- round(haplo.Ord$Total, digits = 8)
haplo.Ord$`data$Dx=1` <- round(haplo.Ord$`data$Dx=1`, digits = 8)
haplo.Ord
colnames(haplo.Ord) <- c('DRB1','DRB4','DQB1', 'DQA1', 'Total_Frequency','Controls_Frequency', 'Cases_Frequency')

# Write EM results 
haplo.Ord['Cases_count'] <- rep(nrow(data.cases), nrow(haplo.Ord))
haplo.Ord['Controls_count'] <- rep(nrow(data.controls), nrow(haplo.Ord))
write.table(haplo.Ord, file = paste(settings$directory$HLA_Output_Haplotype, 'HLA_Haplotype_EM.csv', sep=''),
            quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)

############ COUNT OF HAPLOTYPES OF INTEREST ###########

  # If haplotype type I.
data.haplo1.Hetero <- data %>% filter((HLADRB1_A1 == '07:01' | HLADRB1_A2 == '07:01') &
                               (HLADRB4_A1 == '01:01' | HLADRB4_A2 == '01:01') &
                               (HLADQB1_A1 == '02:02' | HLADQB1_A2 == '02:02') &
                               (HLADQA1_A1 == '02:01' | HLADQA1_A2 == '02:01'))

data.haplo1.Homo <- data %>% filter((HLADRB1_A1 == '07:01' & HLADRB1_A2 == '07:01') &
                                    (HLADRB4_A1 == '01:01' & HLADRB4_A2 == '01:01') &
                                    (HLADQB1_A1 == '02:02' & HLADQB1_A2 == '02:02') &
                                    (HLADQA1_A1 == '02:01' & HLADQA1_A2 == '02:01'))
data.haplo1.Hetero <- data.haplo1.Hetero[which(data.haplo1.Hetero$GWASID %notin% data.haplo1.Homo$GWASID),]

  # If haplotype type II 
data.haplo2.Hetero <- data %>% filter((HLADRB1_A1 == '07:01' | HLADRB1_A2 == '07:01') &
                                 (HLADRB4_A1 == '01:03' | HLADRB4_A2 == '01:03') &
                                 (HLADQB1_A1 == '03:03' | HLADQB1_A2 == '03:03') &
                                 (HLADQA1_A1 == '02:01' | HLADQA1_A2 == '02:01'))

data.haplo2.Homo <- data %>% filter((HLADRB1_A1 == '07:01' & HLADRB1_A2 == '07:01') &
                                      (HLADRB4_A1 == '01:03' & HLADRB4_A2 == '01:03') &
                                      (HLADQB1_A1 == '03:03' & HLADQB1_A2 == '03:03') &
                                      (HLADQA1_A1 == '02:01' & HLADQA1_A2 == '02:01'))
data.haplo2.Hetero <- data.haplo2.Hetero[which(data.haplo2.Hetero$GWASID %notin% data.haplo2.Homo$GWASID),]

  # If mixed haplotype type III
data.intersect <- data %>% filter((HLADRB1_A1 == '07:01' & HLADRB1_A2 == '07:01') &
                                    (HLADRB4_A1 == '01:01' & HLADRB4_A2 == '01:03' | HLADRB4_A1 == '01:03' & HLADRB4_A2 == '01:01') &
                                    (HLADQB1_A1 == '02:02' & HLADQB1_A2 == '03:03' | HLADQB1_A1 == '03:03' & HLADQB1_A2 == '02:02') &
                                    (HLADQA1_A1 == '02:01' & HLADQA1_A2 == '02:01'))
data.haplo1.Hetero <- data.haplo1.Hetero %>% filter(GWASID %notin% data.intersect$GWASID)
data.haplo2.Hetero <- data.haplo2.Hetero %>% filter(GWASID %notin% data.intersect$GWASID)


  # If other 
data.haplo.other <- data %>% filter(GWASID %notin% data.haplo1.Hetero$GWASID & GWASID %notin% data.haplo1.Homo$GWASID &
                                      GWASID %notin% data.haplo2.Hetero$GWASID & GWASID %notin% data.haplo2.Homo$GWASID &
                                      GWASID %notin% data.intersect$GWASID)

# Compute counts
haplo1.Hetero.Count <- data.haplo1.Hetero$Dx %>% table()
haplo1.Homo.Count <- data.haplo1.Homo$Dx %>% table()
haplo2.Hetero.Count <- data.haplo2.Hetero$Dx %>% table()
haplo2.Homo.Count <- data.haplo2.Homo$Dx %>% table(); # None, set to 0 
haplo2.Homo.Count <- c(0,0); names(haplo2.Homo.Count) <- c(0,1)
haplo.intersect.Count <- c(); haplo.intersect.Count['0'] <- 0
haplo.intersect.Count['1'] <- data.intersect$Dx %>% table()
haplo.other.Count <- data.haplo.other$Dx %>% table()

# Compute Frequencies
haplo1.Hetero.Freq <- haplo1.Hetero.Count / c(nrow(data.controls), nrow(data.cases)) * 100
haplo1.Homo.Freq <- haplo1.Homo.Count / c(nrow(data.controls), nrow(data.cases)) * 100
haplo2.Hetero.Freq <- haplo2.Hetero.Count / c(nrow(data.controls), nrow(data.cases)) * 100
haplo2.Homo.Freq <- haplo2.Homo.Count / c(nrow(data.controls), nrow(data.cases)) * 100
haplo.intersect.Freq <- haplo.intersect.Count / c(nrow(data.controls), nrow(data.cases)) * 100
haplo.other.Freq <- haplo.other.Count / c(nrow(data.controls), nrow(data.cases)) * 100

# Compute dataframe 
haploInterest <- data.frame(DRB3_A1 = c('01:01', '01:03', '01:01', '01:01', '01:03', 'X'), 
                            DRB1_A1 = c('07:01', '07:01', '07:01', '07:01', '07:01', 'X'), 
                            DQA1_A1 = c('02:01', '02:01', '02:01', '02:01', '02:01', 'X'), 
                            DQB1_A1 = c('02:02', '03:03', '02:02', '02:02', '03:03', 'X'), 
                            DRB3_A2 = c('01:01', '01:01', '01:03', 'X', 'X', 'X'), 
                            DRB1_A2 = c('07:01', '07:01', '07:01', 'X', 'X', 'X'), 
                            DQA1_A1 = c('02:01', '02:01', '02:01', 'X', 'X', 'X'), 
                            DQB1_A2 = c('02:02', '03:03', '03:03', 'X', 'X','X'), 
                            N_Cases = c(haplo1.Homo.Count['1'], haplo2.Homo.Count['1'],
                                        haplo.intersect.Count['1'], haplo1.Hetero.Count['1'], 
                                        haplo2.Hetero.Count['1'], haplo.other.Count['1']),
                            N_control = c(haplo1.Homo.Count['0'], haplo2.Homo.Count['0'],
                                        haplo.intersect.Count['0'], haplo1.Hetero.Count['0'], 
                                        haplo2.Hetero.Count['0'], haplo.other.Count['0']), 
                            Freq_cases = c(haplo1.Homo.Freq['1'], haplo2.Homo.Freq['1'],
                                        haplo.intersect.Freq['1'], haplo1.Hetero.Freq['1'], 
                                        haplo2.Hetero.Freq['1'], haplo.other.Freq['1']),
                            Freq_contrls = c(haplo1.Homo.Freq['0'], haplo2.Homo.Freq['0'],
                                           haplo.intersect.Freq['0'], haplo1.Hetero.Freq['0'], 
                                           haplo2.Hetero.Freq['0'], haplo.other.Freq['0'])
                            )

# Compute number of cases and controls, initialize loop, get group of reference
ncases <- nrow(data.cases); ncontrols <- nrow(data.controls)
pval <- c(); low_interval <- c(); up_interval <- c(); OR <- c()

# Compute Odd Ratio and Chi squre 
# For each haplotype 
for (i in 1:(nrow(haploInterest)-1)){
  
  # Get data
  haplo <- haploInterest[i,]
  
  # Compute contingincence table, odd ratio and chi sqare
  haplo.table <- matrix(c(haplo$N_Cases, haplo.other.Count['1'], 
                          haplo$N_control, haplo.other.Count['0'] ), nrow = 2)
  oddratio <- oddsratio.wald(haplo.table)
  
  # Append 
  pval <- c(pval, oddratio$p.value['Exposed2','fisher.exact']);
  OR <- c(OR, oddratio$measure['Exposed2','estimate']); 
  low_interval <- c(low_interval, oddratio$measure['Exposed2','lower']); 
  up_interval <- c(up_interval, oddratio$measure['Exposed2','upper'])
}

# Append reference 
pval <- c(pval, 'Reference'); OR <- c(OR, 'Reference');
low_interval <- c(low_interval, 'Reference'); up_interval <- c(up_interval, 'Reference')

# Append to dataframe 
haploInterest$OR <- OR; haploInterest$low_interval <-low_interval; 
haploInterest$up_interval <- up_interval; haploInterest$Fisher.PVAL <- pval

# Write table
write.table(haploInterest, file = paste(settings$directory$HLA_Output_Haplotype, 'HLA_DR7_Haplotype.csv', sep = ''), 
            quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)

############### DR7 specific ##############

# DR7 homozygous, DR7 hetero
dr7.hetero <- data[which(xor(data$HLADRB1_A1 == '07:01', data$HLADRB1_A2 == '07:01')),]
dr7.homo <- data %>% filter(HLADRB1_A1 == '07:01' & HLADRB1_A2 == '07:01')
other.data <- data %>% filter(HLADRB1_A1 != '07:01' & HLADRB1_A2 !=  '07:01')

# Counts and freqs
dr7.homo.cases.counts <- nrow(dr7.homo %>% filter(Dx == 1)); dr7.homo.cases.freq <- dr7.homo.cases.counts / nrow(data.cases) *100
dr7.homo.cntrls.counts <- nrow(dr7.homo %>% filter(Dx == 0)); dr7.homo.cntrls.freq <- dr7.homo.cntrls.counts / nrow(data.controls)*100
dr7.hetero.cases.counts <- nrow(dr7.hetero %>% filter(Dx == 1)); dr7.hetero.cases.freq <- dr7.hetero.cases.counts / nrow(data.cases) *100
dr7.hetero.cntrls.counts <- nrow(dr7.hetero %>% filter(Dx == 0)); dr7.hetero.cntrls.freq <- dr7.hetero.cntrls.counts / nrow(data.controls) *100
other.cases.counts <- nrow(other.data %>% filter (Dx == 1)); other.cases.freqs<- other.cases.counts / nrow(data.cases) *100
other.cntrls.couts <- nrow(other.data %>% filter(Dx == 0)); other.cntrls.freqs <- other.cntrls.couts / nrow(data.controls) *100

# Create dataframe
dr7.haplo <- data.frame(DRB1_A1 = c('07:01', '07:01', 'X'), DRB1_A2 = c('07:01', 'X', 'X'),
                        N_cases = c(dr7.homo.cases.counts, dr7.hetero.cases.counts, other.cases.counts), 
                        N_controls = c(dr7.homo.cntrls.counts, dr7.hetero.cntrls.counts, other.cntrls.couts), 
                        Freq_cases = c(dr7.homo.cases.freq, dr7.hetero.cases.freq, other.cases.freqs), 
                        Freq_controls = c(dr7.homo.cntrls.freq, dr7.hetero.cntrls.freq, other.cntrls.freqs))

# Initialize loop and get reference group 
pval <- c(); low_interval <- c(); up_interval <- c(); OR <- c()

# Compute oddratio
for (i in 1:2){
  
  # Get data 
  haplo <- dr7.haplo[i,]
  
  # Compute contingency table and oddsratio (odds of being DR7 positve and having the disease vs DR7 neg and not having it )
  haplo.table <- matrix(c(haplo$N_cases, other.data %>% filter(Dx == 1) %>% nrow (), 
                          haplo$N_controls, other.data %>% filter(Dx == 0) %>% nrow() ), nrow = 2)
  oddratio <- oddsratio.wald(haplo.table)
  
  # Append 
  pval <- c(pval, oddratio$p.value['Exposed2','fisher.exact']);
  OR <- c(OR, oddratio$measure['Exposed2','estimate']); 
  low_interval <- c(low_interval, oddratio$measure['Exposed2','lower']); 
  up_interval <- c(up_interval, oddratio$measure['Exposed2','upper'])
  
  
}

# Append for group of reference 
OR <- c(OR, 'Reference'); low_interval <- c(low_interval, 'Reference');
up_interval <- c(up_interval, 'Reference'); pval <- c(pval, 'Reference');

# Append to dataframe 
dr7.haplo$OR <- OR; dr7.haplo$low_interval <-low_interval; 
dr7.haplo$up_interval <- up_interval; dr7.haplo$Fisher.PVAL <- pval

# Write 

# Write table
write.table(dr7.haplo, file = paste(settings$directory$HLA_Output_Haplotype, 'HLA_DR7_Count.csv', sep = ''), 
            quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)

