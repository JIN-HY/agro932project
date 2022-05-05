
#ped <- read.table("largedata/newSAP/SAP.chr02.ped", header=F)
#ped$V6 <- 1
#write.table(ped, "largedata/newSAP/SAP.chr02.ped", sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)

fam <- read.table("largedata/newSAP/SAP.chr02.fam", header=F)
fam$V6 <- 1
write.table(fam, "largedata/newSAP/SAP.chr02.fam", sep="\t", row.names=F)


