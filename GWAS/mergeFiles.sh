#!/bin/bash 

# Initialize files and directories
    #DIRS
WD=$(pwd)
DATAGWAS=$(jq -r '.directory.GWAS' GWAS/settingsGWAS.json)
CONTROLSFILTBIN=$(jq -r '.directory.GWAS_controlsFILTBIN' GWAS/settingsGWAS.json)
CASESFILTBIN=$(jq -r '.directory.GWAS_casesFILTBIN' GWAS/settingsGWAS.json)
    #FILES
mergeList=$(jq -r '.file.mergeList' GWAS/settingsGWAS.json)
HU=$(jq -r '.file.HU' GWAS/settingsGWAS.json)

# Go to GWAS directory 
cd $DATAGWAS

# Create temporary folder
mkdir tmp

# Copy controls
cp ${WD}/${CONTROLSFILTBIN}* tmp

# MErge 
cp ${WD}/${CASESFILTBIN}* tmp

# Copy list of files 
cp ${WD}/${mergeList} tmp

# Merge 
cd tmp
bfiles=$(head -n 1 mergeList.txt)
IFS=' ' read -a bfilestr <<< "$bfiles"
IFS='.' read -a bfile <<< "$bfilestr"
awk NR\>1 mergeList.txt > mL.txt
plink --bfile $bfile --merge-list mL.txt --out HU
cd ..

# Move file to GWAS
mv tmp/HU.* $DATAGWAS

# Delete temporary folder
rm -R tmp
