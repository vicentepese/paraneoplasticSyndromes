#!/bin/bash

# Initialize files and directories
WD=$(pwd)
DATAGWAS=$(jq -r '.directory.GWAS' GWAS/settingsGWAS.json)
IbdGenome=$(jq -r '.file.IBD' GWAS/settingsGWAS.json)
cd $DATAGWAS

king -b HU.bed --ibd
mv king.seg ${WD}/$IbdGenome
rm -r *king*
