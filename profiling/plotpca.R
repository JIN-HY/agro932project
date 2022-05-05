pca <- read.table("cache/SAP.chr02.eigenvec", header=T)
pdf("plot/pc1pc2.pdf")
plot(pca$PC1, pca$PC2, xlab="PC1", ylab="PC2")
dev.off()
