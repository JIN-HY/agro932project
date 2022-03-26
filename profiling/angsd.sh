ml anaconda 
conda activate popgen
angsd -vcf-gl largedata/SAP.vcf.gz -nind 10 -fai hg19.fa.gz.fai -domaf 1
