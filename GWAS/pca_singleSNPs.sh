#!/bin/bash

# Go to directory 
cd ../Data

# Count number of subjects 
nsub=$(awk 'END{print NR}' filtHU_singleSNPs.fam)

# Compute PCA
mv modHU_SSNP.fam filtHU_singleSNPs.fam
plink --bfile filtHU_singleSNPs --pca 20 --out HU_singleSNPs

# Move PCs to Resources
mv HU_singleSNPs.eigen* ../GWAS/ResourcesGWAS/
rm HU_singleSNPs.*
