# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)

########## IMPORT ##########
setwd("~/Documents/paraneoplasticSyndromes")

# Import settings
settings <- jsonlite::fromJSON('GWAS/settingsGWAS.json')

 # Read arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  f = settings$file$filtHU_assoc_logistic
} else if (length(args) == 1) {
  if (args[1] == 'singleSNP'){
    f <-settings$file$filtHU_SSNP_assoc_logistcic
  }
}

# Import association file 
assoc.data <- read.table(f, header = TRUE, sep = '', dec = '.')
colnames(assoc.data) <- c("CHR","SNP","BP","A1","TEST","NMISS","OR","STAT","P" )
assoc.data <- assoc.data[which(assoc.data$TEST == 'ADD'),]

######### PLOT ###########

assoc.data <- assoc.data[order(assoc.data$CHR, assoc.data$BP),]

# Data table variables
pvals <- assoc.data$P
x <- c(1:length(pvals))
chr <- assoc.data$CHR
colors <- c()
flag <- 1
for (c in unique(chr)){
  l <- nrow(assoc.data[which(assoc.data$CHR == c),])
  if (flag == 1){
    colors <- c(colors, rep('#F8766D', l))
    flag = 2
  }
  else{
    colors <- c(colors, rep('#00BFC4', l))
    flag = 1
  }
}

# Read cases and control numbers 
filtHU <- read.table('GWAS/Data/filtHU.fam', header = FALSE, sep = '\t')
filtPostMatchList <-read.table(settings$file$filtPostMatchPatList, header = FALSE, sep = ' ')
filtPostMatchList <- filtPostMatchList$V1
filtHU <- filtHU[which(filtHU$V1 %in% filtPostMatchList),]
num.cases <- nrow(filtHU[which(filtHU$V6 == 2),])
num.controls <- nrow(filtHU[which(filtHU$V6 == 1),])

# Create datatable
plot.data <- data.table('pval' =  -log10(pvals), 'pos' =  x, 'colors' = colors, 'chr' = chr)
axisdf <- plot.data %>% group_by(chr) %>% summarize(center=( max(pos) + min(pos)) / 2 )

# Plot
tlt <- paste('GWAS of', settings$diagnosis, 'Number of cases:', as.character(num.cases), 'vs', as.character(num.controls), 'controls', sep = ' ')
if (length(args) == 0){
  png(paste("GWAS/ManhattanPlots/GWAS_",settings$diagnosis,'.png',sep=''), width = 1920, height = 1080)
} else{
  png(paste("GWAS/ManhattanPlots/GWAS_SSNP",settings$diagnosis,'.png',sep=''), width = 1920, height = 1080)
}
ggplot() + geom_point(data = plot.data, aes(pos, pval, color = as.factor(colors)), alpha = 1) +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center ) + 
  labs(x='Chromosome',y='-log10(pval)', title=tlt) +
  geom_hline(yintercept = -log10(5e-8)) + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size=30),
        axis.text  = element_text(size = 20),
        axis.title = element_text(size=25))
dev.off()

