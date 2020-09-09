#!/bin/bash

mv ../Resources/notCmnrsIDs.txt ../Data/GWAS/Controls/
cd ../Data/GWAS/Controls
plink --file HU_controls --no-fid --no-sex --no-parents --exclude notComnrsIDs.txt --make-bed --noweb --out tst 
plink --bfile tst --recode --tab --out pedtst
rm notComnrsIDs.txt

plink --file HU_cases --merge HU_controls.ped HU_controls.map --recode --out merge