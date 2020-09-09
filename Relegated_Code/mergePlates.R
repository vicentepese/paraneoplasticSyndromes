scz <- read.table("~/Downloads/SCZ_MASTER_CASE_CONTROL_parsedHLA_JUL8_2020_2_Digit_Subset_for_Vicente.csv", header = TRUE, sep = ',')
sczOld <-  read.table("~/Downloads/HLA_CALLS_LGL1_CASE_CONTROLS(4).csv", sep = ',', header = TRUE)  # 98 cases
data <- read.table(file = settings$file$HLA_total, header = TRUE, sep = ',')  # 97 cases

sczOld.cases <- sczOld %>% filter(Dx == 1)
GWAS.cases <- sczOld.cases$GWASID

scz$Dx <- rep(0, nrow(scz))
scz$HLADPA1_A1 <- NULL
scz$HLADPA1_A2 <- NULL
data$sample.id <- NULL
cases <- read.table(settings$file$GWASIDsCases) %>% unlist()
controls <- read.table(settings$file$GWASIDsControls) %>% unlist()

# Add diagnosis to data
commonIDs <- intersect(data$GWASID, sczOld$GWASID)
dataMod <- data %>% filter(GWASID %in% commonIDs)
Dx <- c()
for (subj in as.character(dataMod$GWASID)){
  if (subj %in% GWAS.cases){
    Dx <- c(Dx, 1)
  } else{
    Dx <- c(Dx, 0)
  } 
}
dataMod$Dx <- Dx
dataMod$sample.id <- NULL


finalPlates <- rbind(dataMod, scz)
finalPlates <- rbind(finalPlates, sczOld.cases)

table(finalPlates$Dx)

write.table(finalPlates, settings$file$HLA_calls_all_plates, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
