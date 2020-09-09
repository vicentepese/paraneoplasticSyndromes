#!/bin/bash

########### CONTROLS #############
cd Data/GWAS/Controls
controls="HU_controls" 

# Convert to binary 
files=$(ls)
if [[ $files == *"HU_controls.ped"* ]] 
then 
    if [[ $files != *"BHU_controls.bed"* ]]
    then 
        echo "Creating duplicate vars"
        plink --file $controls --no-sex --no-parents --no-fid --list-duplicate-vars --noweb --out BHU_controls
        echo "Converting HU_controls.ped to binary"
        plink --file $controls --no-fid --no-sex --no-parents --noweb --make-bed --out BHU_controls
        mv BHU_controls.log Log/BHU_controls.log
        echo "Conversion to binary was successful"
        plink --bfile 
    else
        echo "HU_controls already convervet to binary"
    fi 
else
    echo "HU_controls.ped not found. Please run getControls.py to create the .ped file"
fi 

# Remove duplicates
echo "Removing duplicates of HU_controls"
plink --file HU_controls --no-sex --no-parents --no-fid --list-duplicate-vars suppress-first --noweb --out BHU_controls
mv BHU_controls.dupvar Log/BHU_controls.dupvar
echo "Duplicates removed"



# Perform quality controls 
files=$(ls)
log=$(ls)
if [[ $files == *"BHU_controls.bed"* ]] && [[ $log =~ *"BHU_controls.fmendel"* ]]
then 
    echo "Performing quality control"
    echo "Missing genotypes"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --missing --out BHU_controls
    echo "Missingness of case/control status"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --test-missing --out BHU_controls
    echo "Haplotype-based test for non-random missing genotype data"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --test-mishap  --out BHU_controls
    echo "Hardy-Weinberg Equilibrium"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --hardy  --out BHU_controls
    echo "Allele frequency"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --freq  --out BHU_controls
    echo "Mendel error"
    plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --mendel  --out BHU_controls

else
    echo "File not converted to binary. Please convert to binary before performing quality control."
fi 

########### CASES #############
cd ~/Documents/HU
cd Data/GWAS/Cases
cases="HU_cases" 

# Convert to binary 
files=$(ls)
if [[ $files == *"HU_cases.ped"* ]] 
then 
    if [[ $files != *"BHU_cases.bed"* ]]
    then 
        echo "Converting HU_cases.ped to binary"
        plink --file $cases --no-fid --no-sex --no-parents --noweb --make-bed --out BHU_cases
        echo ""
        mv BHU_cases.log Log/BHU_cases.log
        echo "Conversion to binary was successful"
    else
        echo "HU_cases already convervet to binary"
    fi 
else
    echo "HU_cases.ped not found. Please run getCases.py to create the .ped file"
fi 

# Perform quality controls 
files=$(ls)
log=$(ls)
if [[ $files == *"BHU_cases.bed"* ]]
then 
    echo "Performing quality control"
    echo "Missing genotypes"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --missing --out Log/BHU_cases
    echo "Missingness of case/control status"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --test-missing --out Log/BHU_cases
    echo "Haplotype-based test for non-random missing genotype data"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --test-mishap  --out Log/BHU_cases
    echo "Hardy-Weinberg Equilibrium"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --hardy  --out Log/BHU_cases
    echo "Allele frequency"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --freq  --out Log/BHU_cases
    echo "Mendel error"
    plink --bfile BHU_cases --no-fid --no-sex --no-parents --noweb --mendel  --out Log/BHU_cases

else
    echo "File not converted to binary. Please convert to binary before performing quality control."
fi 

# plink --file $controls --no-fid --no-sex --no-parents --noweb --make-bed --out BHU_controls

# Missing genotypes
#plink --bfile BHU_controls --no-fid --no-sex --no-parents --noweb --missing

# Allele frequency statistics
# plink --bfile BForward --freq --out freq_stat --noweb

#Sort MAF
# sort --key=5 -nr freq_stat.frq > sortFreq_stat.frq
# awk '$5 < 0.05 && $5 > 0 {print $0}' sortFreq_stat.frq > filtFreq_stat.frq 

# Missing statistics
# plink --bfile BForward --missing --out miss_stat --noweb

#  Association analysis 
# plink --bfile BForward --assoc --out ../../Output/ForwardAssoc --noweb

plink --file BHU_cases --no-sex --no-fid --no-parents --merge HU_controls.ped HU_controls.map --make-bed tst 
plink --file HU_cases --no-sex --no-fid --no-parents --make-bed --out BHU_cases
plink --file HU_controls --no-sex --no-fid --no-parents --make-bed --out BHU_controls
plink --bfile BHU_cases --bmerge BHU_controls.bed BHU_controls.bim BHU_controls.fam --make-bed --out merge