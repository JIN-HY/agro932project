ml anaconda 
conda activate popgen
angsd -vcf-gl largedata/SAP.vcf.gz -nind 358 -fai largedata/Sbicolor_454_v3.0.1.fa.fai -out cache/output -doMajorMinor 1 -doMaf 8 -doCounts 0 -doSaf 2 -uniqueOnly 0 -anc largedata/Sbicolor_454_v3.0.1.fa -ref largedata/Sbicolor_454_v3.0.1.fa
#angsd -bam bamlist.txt -out output -doMajorMinor 1 -doMaf 1 -doSaf 2 -uniqueOnly 0 -anc Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -minMapQ 30 -minQ 20 -nInd 20 -baq 1 -ref Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -GL 1
# use realSFS to calculate sfs
#realSFS output.saf.idx -fold 1 > output.sfs
