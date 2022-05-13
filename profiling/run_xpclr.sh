#!/bin/bash
ml xpclr
ml tabix

xpclr --input SAP.chr02.vcf.gz --out ../../cache/xpclr_res.txt --format vcf --samplesA ../../cache/chr2grp1.txt --samplesB ../../cache/chr2grp2.txt --ld 0.7 --size 100000 --step 10000 --chr Chr02
