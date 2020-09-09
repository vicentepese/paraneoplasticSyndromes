#!/bin/bash 

# Initialize directories and files 
    # DIRS
WD=$(pwd)
DATAGWAS=$(jq -r '.directory.GWAS' GWAS/settingsGWAS.json)
DATA=$(jq -r '.directory.Data' GWAS/settingsGWAS.json)
    # FILES
excludeID=$(jq -r '.file.excludeID' GWAS/settingsGWAS.json)
HU=$(jq -r '.file.HU' GWAS/settingsGWAS.json)
FILTHU=$(jq -r '.file.filtHU' GWAS/settingsGWAS.json)


# Go to directory 


# Rewrite files with 10% missingness of genotype and sample, and remove MHC region in CHR6 
plink --bfile $WD/$FILTHU --remove $WD/$excludeID\
    --no-fid --no-sex --no-parents --not-chr 25,26 \
    --maf 0.05 --geno 0.1 --mind 0.1 \
    --make-bed --out $WD/$FILTHU
