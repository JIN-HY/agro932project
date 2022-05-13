#!/bin/bash

ml vcftools

#vcftools --bcf snps.bcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --fst-window-size 10000 --fst-window-step 1000 --out win_1k
vcftools --gzvcf largedata/newSAP/SAP.chr02.vcf.gz --out cache/chr2tajima --TajimaD 100000
