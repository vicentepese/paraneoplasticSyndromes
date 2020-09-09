#!/bin/bash 

# Initialize files and directories 
    #DIRS
WD=$(pwd)
DATAGWAS=$(jq -r '.directory.GWAS' GWAS/settingsGWAS.json)
DATA=$(jq -r '.directory.Data' GWAS/settingsGWAS.json)
    #FILES
excludeID=$(jq -r '.file.excludeID' GWAS/settingsGWAS.json)
HU=$(jq -r '.file.HU' GWAS/settingsGWAS.json)
FILTHU=$(jq -r '.file.filtHU' GWAS/settingsGWAS.json)

# Go to directory 
cd $DATAGWAS

# Rewrite files with 10% missingness of genotype and sample, and remove MHC region in CHR6 
plink --bfile $WD/$HU --no-fid --no-sex --no-parents --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp
rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
rm -r temp*
plink --bfile $WD/$HU --remove ${WD}/$excludeID  \
    --no-fid --no-sex --no-parents \
    --exclude-snps $rsExclude --not-chr XY \
    --maf 0.05 --geno 0.1 --mind 0.15 \
    --make-bed --out $WD/$FILTHU

mv filtHU* $DATA
