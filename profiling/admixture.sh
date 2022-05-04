#!/bin/bash

ml anaconda
conda activate popgen
#for K in {1..8};
#do \
admixture --cv largedata/SAPprune.bed $1 | tee log$1.out;
#done
#grep -h CV log*.out
#admixture largedata/SAP.bed 6
