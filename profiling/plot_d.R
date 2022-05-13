tajimad <- read.table("cache/chr2tajima.Tajima.D", header=T)
pdf("plot/tajimad.pdf")
plot(tajimad$BIN_START, tajimad$TajimaD, type="l", xlab="Posistion (bp)", ylab="Tajima's D")
dev.off()
