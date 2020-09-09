#!/bin/bash 

# Go to directory 
cd ../Data/

# Rewrite files with 10% missingness of genotype and sample, and remove MHC region in CHR6 
plink --bfile filtHU --no-fid --no-sex --no-parents --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out temp
rsExclude=$(awk '{ ORS=", "};{print $2}' temp.bim)
rm -r temp*
plink --bfile filtHU \
    --no-fid --no-sex --no-parents \
    --exclude-snps $rsExclude --not-chr XY \
    --maf 0.05 --geno 0.1 --mind 0.1 \
    --make-bed --out filtHU_singleSNPs

