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
