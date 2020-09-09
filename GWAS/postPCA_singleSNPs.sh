#!/bin/bash 

cd ../Data
plink --bfile filtHU_singleSNPs --keep ../GWAS/ResourcesGWAS/patList_SSNP.txt --make-bed --out temp
plink --bfile temp --pca 20 --out postHU_singleSNPs
mv postHU_singleSNPs.eigen* ../GWAS/ResourcesGWAS/
rm -r postHU_singleSNPs*
rm -r temp*