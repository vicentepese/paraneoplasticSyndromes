#!/bin/bash 

cd ../Data
plink --bfile filtHU --keep ../GWAS/ResourcesGWAS/patList.txt --make-bed --out temp
plink --bfile temp --pca 20 --out postHU
mv postHU.eigen* ../GWAS/ResourcesGWAS/
rm -r postHU*
rm -r temp*