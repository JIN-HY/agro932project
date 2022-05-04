#!/bin/bash
ml plink
plink --bfile largedata/SAP --indep-pairwise 50 10 0.1
plink --bfile largedata/SAP --extract plink.prune.in --make-bed --out largedata/SAPprune
