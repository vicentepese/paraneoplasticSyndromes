#!/bin/bash

# Initialize directories and files
    #DIRS
DATA=$(jq -r '.directory.Data' GWAS/settingsGWAS.json)
RESOURCESGWAS=$(jq -r '.directory.ResourcesGWAS' GWAS/settingsGWAS.json)
WD=$(pwd)
FILTHU=$(jq -r '.file.filtHU' GWAS/settingsGWAS.json)
    # FILES
modHU_FAM=$(jq -r '.file.modHU' GWAS/settingsGWAS.json)
filtHU_FAM=$(jq -r '.file.filtHU_FAM' GWAS/settingsGWAS.json)

# Go to directory  
# Count number of subjects 
nsub=$(awk 'END{print NR}' $filtHU_FAM)

# Compute PCA
mv $modHU_FAM $filtHU_FAM
plink --bfile $FILTHU --pca 20 --out HU

# Move PCs to Resources
mv HU.eigen* ${WD}/$RESOURCESGWAS
rm HU.*
