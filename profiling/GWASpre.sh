
ml plink vcftools gemma

plink --vcf SAP.chr02.vcf.gz --make-bed --out SAP.chr02
plink -bfile SAP.chr02 --freq --missing --out SAP.chr02

#vcftools --gzvcf SAP.chr02.vcf.gz --max-missing 0.2 --maf 0.04 --recode --out SAP.chr02
#gzip SAP.chr02.recode.vcf
#plink --vcf SAP.chr02.recode.vcf.gz --make-bed --out SAP.chr02.recode
plink -bfile SAP.chr02 --pca 'header' --out SAP.chr02

gemma -bfile SAP.chr02 -gk 1 -o SAP.chr02

gemma -bfile SAP.chr02 -c pc3.txt -k output/SAP.chr02.recode.cXX.txt -p phenoHN.txt -lmm 4 -n ??? -o FloweringTimeHN -miss 0.8 -r2 1 -hwe 0 -maf 0.04
