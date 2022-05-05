result <- read.table("cache/SAPprune.5.Q")
result$grp <- max.col(result)
result.grp <- list()
for (g in unique(result$grp)){
  result.grp[[g]] <- result[result$grp==g,]
  result.grp[[g]] <- result.grp[[g]][order(result.grp[[g]][,g], decreasing = T),]
}
result <- do.call("rbind", result.grp)
result1 <- result[, -6]
#result1 <- result[with(result,order(grp,V1,V2,V3,V4,V5,decreasing = F)),-6]
# STILL NEED TO FIX ORDER DONT NEED ANY MORE
pdf("plot/admixture5.pdf", height = 3)
barplot(t(as.matrix(result1)), col=rainbow(5), xlab="Individual", ylab="Ancestray", border=NA, names.arg = c(1:358))
dev.off()

genotypes <- read.table("largedata/SAP.fam")
g_result <- merge(genotypes, result, by = 0)
g_result[g_result$V1.x == "PI_564163", ]
g_result$bicolor <- g_result$grp==4

write.table(g_result, "cache/subpopulation.csv", quote = F, sep = ",", row.names=F)


pca <- read.table("cache/SAP.chr02.eigenvec", header=T)
pdf("plot/pc1pc2.pdf")
plot(pca$PC1, pca$PC2, xlab="PC1", ylab="PC2")
outliers <- pca[pca$PC1 < -0.2,]
text(outliers$PC1, outliers$PC2, labels = outliers$FID, adj = c(0,0))
dev.off()

outliergrp <- g_result[gsub('_', '', g_result$V1.x) %in% outliers$FID,]
write.table(outliergrp[,c(2,8:13)], "cache/outliers.csv", quote = F, sep = ",", row.names=F)
