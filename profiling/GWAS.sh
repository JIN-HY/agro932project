
ml plink vcftools gemma R
cd largedata/newSAP
plink --vcf SAP.chr02.vcf.gz --make-bed --out SAP.chr02
plink -bfile SAP.chr02 --freq --missing --out SAP.chr02

#vcftools --gzvcf SAP.chr02.vcf.gz --max-missing 0.2 --maf 0.04 --recode --out SAP.chr02
#gzip SAP.chr02.recode.vcf
#plink --vcf SAP.chr02.recode.vcf.gz --make-bed --out SAP.chr02.recode
plink -bfile SAP.chr02 --pca 'header' --out SAP.chr02
cp SAP.chr02.eigenvec ../../cache
cd ../../cache
cut -d' ' -f3-5 SAP.chr02.eigenvec | tail -n +2 | awk '{print 1, $0}' > pc3.txt
sed -i 's/ /\t/g' pc3.txt

Rscript profiling/correct-9.R
gemma -bfile SAP.chr02 -gk 1 -o SAP.chr02


# phenotypic 
cd ..
sed 's/,/\t/g' data/DaysToBloom_HNblup.csv | sed 's/_//g' | sed 's/"//g' |tail -n +2> cache/phenoHN.txt
sed 's/,/\t/g' data/DaysToBloom_LNblup.csv | sed 's/_//g' | sed 's/"//g' | tail -n +2> cache/phenoLN.txt

# GWAS in home directory
gemma -bfile largedata/newSAP/SAP.chr02 -c cache/pc3.txt -k cache/SAP.chr02.cXX.txt -p cache/phenoHN.txt -lmm 4 -n 2 -o FloweringTimeHN -miss 0.1 -r2 1 -hwe 0 -maf 0.04
gemma -bfile largedata/newSAP/SAP.chr02 -c cache/pc3.txt -k cache/SAP.chr02.cXX.txt -p cache/phenoLN.txt -lmm 4 -n 2 -o FloweringTimeLN -miss 0.1 -r2 1 -hwe 0 -maf 0.04
