#!/bin/bash 

# Compute distance between cases 
cd ../Data 
plink --bfile filtHU --distance square --out distance
mv *distance* ../GWAS/ResourcesGWAS/