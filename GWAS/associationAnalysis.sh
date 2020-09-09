#!/bin/bash

# Initialize files and directories
    #DIRS
    DATA=$(jq -r '.directory.Data' GWAS/settingsGWAS.json)
    WD=$(pwd)
    #FILES
    filtPostMatchPatList=$(jq -r '.file.filtPostMatchPatList' GWAS/settingsGWAS.json)
    PCS=$(jq -r '.file.PCS' GWAS/settingsGWAS.json)
    filtHU_assoc_logistic=$(jq -r '.file.filtHU_assoc_logistic' GWAS/settingsGWAS.json)
    modHUFAM=$(jq -r '.file.modHU_fam' GWAS/settingsGWAS.json)
    filtHUFAM=$(jq -r '.file.filtHU_fam' GWAS/settingsGWAS.json)
    filtHU=$(jq -r '.file.filtHU' GWAS/settingsGWAS.json)

# Go to directory 
cd $DATA

# Copy genotyped fam file 

cp $WD/$modHU_fam $WD/$filtHUFAM
plink --bfile $WD/$filtHU --keep ${WD}/$filtPostMatchPatList \
    --covar ${WD}/$PCS --covar-name PC1, PC2, PC3, PC4 \
    --logistic --allow-no-sex --out $WD/$filtHU
mv $WD/${filtHU}.assoc.logistic ${WD}/$filtHU_assoc_logistic
