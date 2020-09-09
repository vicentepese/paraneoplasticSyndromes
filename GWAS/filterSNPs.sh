#!/bin/bash

# Initialize directories and files 
    #DIRS
WD=$(pwd)
GWAS_controlsBIN=$(jq -r '.directory.GWAS_controlsBIN' GWAS/settingsGWAS.json)
GWAS_casesBIN=$(jq -r '.directory.GWAS_casesBIN' GWAS/settingsGWAS.json)
GWAS_controlsFILTBIN=$(jq -r '.directory.GWAS_controlsFILTBIN' GWAS/settingsGWAS.json)
GWAS_casesFILTBIN=$(jq -r '.directory.GWAS_casesFILTBIN' GWAS/settingsGWAS.json)
    #FILES
commonSNPs=$(jq -r '.file.commonSNPs' GWAS/settingsGWAS.json)


########## CONTROLS ############

# Copy directory 
cp $(echo $commonSNPs) "${GWAS_controlsBIN}"commonSNPs.txt

# Go to binaries directory
cd $GWAS_controlsBIN

files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bed"* ]] ; then 
        echo "Filtering SNPs in " $file
        IFS='.' read -a strarr <<< "$file"
        plink --bfile ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --extract commonSNPs.txt --make-bed --out $WD/$GWAS_controlsFILTBIN${strarr[0]}
        mv ${strarr[0]}.log $WD/${GWAS_casesFILTBIN}/Log/${strarr[0]}.log
    fi 
done
rm commonSNPs.txt

cd $GWAS_controlsFILTBIN
files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bim"* ]] ; then 
        echo $file
        awk 'END{print NR}' $file
    fi 
done 

########### CASES ##############

cd $WD
cp $commonSNPs "$GWAS_casesBIN"/commonSNPs.txt

# Go to binaries directory
cd $GWAS_casesBIN

files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bed"* ]] ; then 
        echo "Filtering SNPs in " $file
        IFS='.' read -a strarr <<< "$file"
        plink --bfile ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --extract commonSNPs.txt --make-bed --out $WD/${GWAS_casesFILTBIN}/${strarr[0]}
        mv ${strarr[0]}.log $WD/${GWAS_casesFILTBIN}/Log/${strarr[0]}.log
    fi 
done
rm commonSNPs.txt

cd $GWAS_casesFILTBIN
files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bim"* ]] ; then 
        echo $file
        awk 'END{print NR}' $file
    fi 
done 