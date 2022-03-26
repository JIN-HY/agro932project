#!/bin/bash

ml tassel samtools
run_pipeline.pl -h largedata/SAP_imputed.hmp -export largedata/SAP.vcf -exportType VCF
bgzip largedata/SAP.vcf.gz
tabix -p vcf largedata/SAP.vcf.gz

