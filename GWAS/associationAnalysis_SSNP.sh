# Read directories / files
DATADIR=`jq '.directory.Data' settingsGWAS.json` 
filtHU_SSNP_assoc_logistic=`jq '.file.filtHU_SSNP_assoc_logistcic' settingsGWAS.json`


# Go to Data directory
cd $DATADIR

# Compute association analysis
plink --bfile filtHU_singleSNPs --keep ResourcesGWAS/filtPostMatchPatList_SSNP.txt \
    --covar ResourcesGWAS/PCAtst_SSNP.eigenvec --covar-name PC1, PC2, PC3, PC4 \
    --logistic --allow-no-sex --out filtHU_SSNP
mv filtHU_SSNP.assoc.logistic $filtHU_SSNP_assoc_logistic
