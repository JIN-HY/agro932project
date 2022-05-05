
pheHN <- read.csv("cache/phenoHN.txt", sep = "\t", header = F)
pheLN <- read.csv("cache/phenoLN.txt", sep = "\t", header = F)

fam <- read.csv("largedata/newSAP/SAP.chr02.fam", sep = "\t", header = F)
#fam <- fam[,1]

pheHN <- merge(pheHN, fam, all.y = T, by.x = "V1", by.y = "V1")[,c(1,2)]
pheLN <- merge(pheLN, fam, all.y = T, by.x = "V1", by.y = "V1")[,c(1,2)]

pheHN[is.na(pheHN[,2]),2]="NA"
pheLN[is.na(pheHN[,2]),2]="NA" 

pheHN <- pheHN[match(fam[,1], pheHN[,1]),]
pheHN <- pheHN[match(fam[,1], pheHN[,1]),]

write.table(pheHN, "cache/phenoHN.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(pheLN, "cache/phenoLN.txt", sep = "\t", quote = F, row.names = F, col.names = F)

