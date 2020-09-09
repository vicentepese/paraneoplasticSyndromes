#!/bin/bash 

########## CONTROLS ###########

# - No missalignment in controls. No "copys". This section takes each .ped and .map file and 
#       converst it to a binary that is then moved to the Binaries folder
# - Logs are moved to Logs directory

# Get files and direcotories
    #DIRS
WD=$(pwd)
CONTROLSOG=$(jq -r '.directory.GWAS_controlsOG' GWAS/settingsGWAS.json)
CONTROLSBIN=$(jq -r '.directory.GWAS_controlsBIN' GWAS/settingsGWAS.json)

# Go to controls 
cd $CONTROLSOG

# Convert controls to binarie
files=($(ls))
for file in ${files[@]} ; do
    if [[ $file == *'.ped'* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents \
            --noweb --make-bed --out $WD/${CONTROLSBIN}${strarr[0]} >> $WD/${CONTROLSBIN}${strarr[0]}.log
        mv $WD/${CONTROLSBIN}${strarr[0]}.log $WD/${CONTROLSBIN}Log/${strarr[0]}.log
    fi  
done 
cd $WD

########## CASES ###########

# There is are copys of cases in the .ped file. 
# - .ped files are filtered and copies are rejected. 
# - Files are converted to binaries, sent to the ../Binaries directory and .log to the ../Log directory

# Get files and directories
    # DIRS
CASESOG=$(jq -r '.directory.GWAS_casesOG' GWAS/settingsGWAS.json)
CASESBIN=$(jq -r '.directory.GWAS_casesBIN' GWAS/settingsGWAS.json)

# Go to cases
cd $CASESOG

# Filter data
files=($(ls))
if [[ $files != *'_filt.ped'* ]] ; then 
    files=($(ls))
    for file in ${files[@]} ; do 
        if [[ $file == *'.ped'* ]] ; then 
            echo "Removing duplicate patients of " $file 
            IFS='.' read -a strarr <<< "$file"  
            awk '$3 !~ /Copy/ { print $0 ; }' $file > ${strarr[0]}"_filt.ped"
            cp ${strarr[0]}.map ${strarr[0]}_filt.map
            echo "Duplicates removed"
        fi 
    done
else 
    echo "Cases already filtered"
fi 


files=($(ls))
for file in ${files[@]} ; do
    if [[ $file == *'.ped'* ]] && [[ $file == *"filt"* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --noweb --make-bed --out $WD/$CASESBIN${strarr[0]} >> $WD/$CASESBIN${strarr[0]}.log
        mv $WD/$CASESBIN/${strarr[0]}.log $WD/${CASESBIN}Log/${strarr[0]}.log
    fi 
done 

# Remove
rm -R *filt*