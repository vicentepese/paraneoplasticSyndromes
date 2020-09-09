# Import libraries
library(xlsx)
library(readxl)
library(jsonlite)
library(dplyr)
library(stringr)

############ INITIALIZATION ############ 

# Set working directory 
setwd("~/Documents/paraneoplasticSyndromes/")

# Import settings
settings <- jsonlite::read_json('settings.json')

# Get args
args = commandArgs(trailingOnly=TRUE)
if (length(args) ==1){
  settings$diagnosis <- args[1]
}

# Import data
paraneo.data <- read_xlsx(settings$file$GWAS_total)
table(paraneo.data$Dx)
unique(paraneo.data$Dx)
table(paraneo.data$`Dx-ling`)
unique(paraneo.data$`Dx-ling`)

########### GET PLATES CONTROLS ###############
print('Parsing controls')
# Get controls data and count 
ctrls.dx <- c("NMDA control", "HIMC", "F or G3", "G3", "G3 or J", "sleep study", "k", "APOE study", "control", "K", 
  "Control", "Healthy Control", "Normal", "PSG Sstudy", "PSG study", "normal")
ctrls.data <- paraneo.data[which(paraneo.data$Dx %in% ctrls.dx),]
table(ctrls.data$Dx)

# Get plates
filterplates <- function(GWASID){
  return(strsplit(gsub("([0-9]*)([A-Z]*)([0-9]*)", "\\1 \\2 \\3", GWASID), " ")[[1]][1])
}

ctrls.plates <- sapply(ctrls.data$`GWAS ID`, FUN = filterplates)
table(ctrls.plates)
ctrls.plates <- as.numeric(unique(ctrls.plates))
ctrls.plates <- ctrls.plates[which(ctrls.plates >= 77)]
ctrls.plates

# Get control files from second batch (plates 120 and 121)
kls.data <- read_xlsx(settings$file$GWAS_KLS)
kls.GWASID <- kls.data[which(kls.data$Dx %in% ctrls.dx),]$`GWAS ID` 

# Write 
write.table(as.numeric(ctrls.plates), file = settings$file$plates, row.names = FALSE, 
          col.names = FALSE, sep = ',')

# Get GWAS ID
GWAS.ID <- c(ctrls.data$`GWAS ID`, kls.GWASID) %>% unique()
GWAS.ID <- GWAS.ID[which(!is.na(GWAS.ID))]

# Write
write.table(GWAS.ID, file = settings$file$GWASIDsControls, row.names =  FALSE,
            col.names = FALSE, sep = ',')

print('Controls parsed')
print(paste('Control GWAS Ids saved in', settings$file$GWASIDsControls, sep = ' '))

########### GET PLATES CASES ##########

print('Parsing cases')
print(paste('Diagnosis', settings$diagnosis, sep = ' '))

# Get cases data 
switch(settings$diagnosis,
       
       'HU'={
         paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
         table(paraneo.cases$Dx)
         cases.data <- paraneo.cases[which(paraneo.cases$Dx == 'HU'),]
         
         # Get GWAS ID
         GWAS.ID.cases <- cases.data$`GWAS ID` %>% unique()
         
         # Write 
         write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
                     col.names = FALSE, sep = ',')
       },
       
       'Caspr'={
         # Initialize
         GWAS.ID.cases <- c()
         
         # From paraneo.cases
         paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
         table(paraneo.cases$Dx)
         GWAS.ID.cases <- paraneo.cases[which(paraneo.cases$Dx == 'Caspr2' | paraneo.cases$Dx == 'CASPR2'),]$`GWAS ID` 
         
         # From paraneo.data
         GWAS.ID.cases <- c(GWAS.ID.cases, paraneo.data[which(paraneo.data$Dx == 'Caspr2' | paraneo.data$Dx == 'CASPR2'),]$`GWAS ID` )
         
         # Unique 
         GWAS.ID.cases <- GWAS.ID.cases %>% unique()
         
         # Write 
         write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
                     col.names = FALSE, sep = ',')
         
       },
       
       'LGI1'={
         paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
         table(paraneo.cases$Dx)
         cases.data <- paraneo.cases[which(paraneo.cases$Dx == 'LGI1'),]
         
         # Get plates 
         cases.plates <- sapply(cases.data$`GWAS ID`, FUN = filterplates)
         table(cases.plates)
         unique(cases.plates)
         
         
         # Get GWAS ID
         GWAS.ID.cases <- cases.data$`GWAS ID` %>% unique()
         
         # Write 
         write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
                     col.names = FALSE, sep = ',')
       },
       
       'YO'={
         paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
         cases.data <- paraneo.cases[which(paraneo.cases$Dx == 'Yo' | paraneo.cases$Dx == 'YO'),]
         
         # Get GWAS ID
         GWAS.ID.cases <- cases.data$`GWAS ID` %>% unique()
         
         # Write 
         write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
                     col.names = FALSE, sep = ',')
       },
       
       'AK5'={
         paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
         cases.data <- paraneo.cases[which(paraneo.cases$Dx == 'AK5'),]
         
         # Get GWAS ID
         GWAS.ID.cases <- cases.data$`GWAS ID` %>% unique()
         
         # Write 
         write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
                     col.names = FALSE, sep = ',')
       })

print('Cases parsed')
print(paste('Cases GWAS Ids saved in', settings$file$GWASIDsCases, sep  = ' '))
print('Patient information successfully parsed')

# Write tables
cases.table <- table(paraneo.cases$Dx)
write.table(cases.table, settings$file$Cases_samples, sep = ',', row.names = FALSE, col.names = c('Disease', 'Number of cases'), quote = FALSE)
controls.table <- table(paraneo.data$Dx)
write.table(controls.table, settings$file$Controls_samples, sep = ',', row.names = FALSE, col.names = c('Disease', 'Number of cases'), quote = FALSE)

