
ml vcftools plink
vcftools --gzvcf largedata/SAP.vcf.gz --plink --out largedata/SAP
plink --file largedata/SAP --make-bed --out largedata/SAP
