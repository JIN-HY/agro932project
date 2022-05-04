
library(MASS)
library(candisc)
library(e1071)
library(g3tools)
library(lme4)

SAP <- read.csv("SAPstat873.csv")
write.table(SAP, "SAPstat873.csv",sep = ",",quote = F,row.names = F)
#SAP$`PoorStand.`=as.factor(SAP$`PoorStand.`)
SAP$`PoorStand.`=NULL
SAP$Treatment <- factor(SAP$Treatment)
SAP$Block <- factor(SAP$Block)

SAP.HN <- SAP[SAP$Treatment=='SufficientNitrogen',]
SAP.LN <- SAP[SAP$Treatment=='LowNitrogen',]

H2HN <- data.frame(trait=colnames(SAP)[9:ncol(SAP)], H2=0)
H2LN <- data.frame(trait=colnames(SAP)[9:ncol(SAP)], H2=0)

for (t in 9:ncol(SAP.HN)) {
  trait=colnames(SAP.HN)[t]
  fit=get_BLUP(data=SAP.HN,model = as.formula(paste(trait,"~(1|SorghumAccessionRed) + (1|Block)")),which.factor="SorghumAccessionRed",outfile=paste0(trait,"_HNblup.csv"))
  #print(summary(fit))
  H2HN[t-2,2]=(get_H2(fit,numerator = "GenotypeID",denominator = data.frame(f=c("GenotypeID","Residual"),df=c(1,2))))
}

for (t in 9:ncol(SAP.LN)) {
  trait=colnames(SAP.LN)[t]
  fit=get_BLUP(data=SAP.LN,model = as.formula(paste(trait,"~(1|SorghumAccessionRed) + (1|Block)")),which.factor="SorghumAccessionRed",outfile=paste0(trait,"_LNblup.csv"))
  #print(summary(fit))
  H2LN[t-2,2]=(get_H2(fit,numerator = "GenotypeID",denominator = data.frame(f=c("GenotypeID","Residual"),df=c(1,2))))
}

# for (t in 9:ncol(SAP)) {
#   trait=colnames(SAP)[t]
#   fit=get_BLUP(data=SAP,model = as.formula(paste(trait,"~(1|SorghumAccessionRed) + Treatment")),which.factor="SorghumAccessionRed",outfile=paste0(trait,"_blup.csv"))
#   #print(summary(fit))
#   #H2s[t-2,2]=(get_H2(fit,numerator = "GenotypeID",denominator = data.frame(f=c("GenotypeID","Residual"),df=c(1,2))))
# }

blupLN=data.frame(SorghumAccessionRed = unique(SAP$SorghumAccessionRed), Treatment = 'LN')
blupHN=data.frame(SorghumAccessionRed = unique(SAP$SorghumAccessionRed), Treatment = 'HN')
#colnames(blupdfs)="Taxa"
for (t in 9:ncol(SAP)) {
  trait=colnames(SAP)[t]
  blupfile=read.csv(paste0(trait,"_HNblup.csv"))
  blupHN=merge(blupHN,blupfile,by.x = "SorghumAccessionRed",by.y="X",all.x = T)
  colnames(blupHN)[t-6]=trait
}
for (t in 9:ncol(SAP)) {
  trait=colnames(SAP)[t]
  blupfile=read.csv(paste0(trait,"_LNblup.csv"))
  blupLN=merge(blupLN,blupfile,by.x = "SorghumAccessionRed",by.y="X",all.x = T)
  colnames(blupLN)[t-6]=trait
}

SAPblup = rbind(blupLN, blupHN)
SAPblup = na.omit(SAPblup)
write.table(SAPblup,"SAPnitrogenblup2020.csv",sep = ",",quote = F,row.names = F)



for (col in 9:26) {
  print(col)
  qqnorm(SAP[,col])
  sk <- skewness(SAP[,col], na.rm = T)
  if (sk >0.8 & sk<1.2){
    sk <- skewness(log(SAP[,col]), na.rm = T)
  } else if (sk >1.2 & sk <3){
    sk <- skewness(sqrt(SAP[,col]), na.rm =T)
  } else if (sk>=3){
    sk <- skewness(SAP[,col]^(1/3), na.rm = T)
  }
  print(sk)
}

SAP.trans <- SAPblup
for (col in 3:20) {
  sk <- skewness(SAPblup[,col], na.rm = T)
  if (sk >0.8 & sk<1.2){
    SAP.trans[,col] <- log(SAP.trans[,col])
  } else if (sk >1.2 & sk<3){
    SAP.trans[,col] <- sqrt(SAPblup[,col])
  } else if (sk>=3){
    SAP.trans[,col] <- SAPblup[,col]^(1/3)
  }
}


SAP.mod = lm(cbind(DaysToBloom,MedianLeafAngle,LeafAngleSDV,PaniclesPerPlot,
                     PanicleGrainWeight,EstimatedPlotYield,FlagLeafLength,
                     FlagLeafWidth,ExtantLeafNumber,PlantHeight,TillersPerPlant,
                     StemDiameterLower,StemDiameterUpper,RachisLength,
                     RachisDiameterLower,RachisDiameterUpper,PrimaryBranchNo,BranchInternodeLength)
               ~Treatment+SorghumAccessionRed,data=SAP.trans)

SAP.anova.out=Anova(SAP.mod,type='III',test="Wilks")
SAP.can.t=candisc(SAP.mod,term = "Treatment")
SAP.can.g=candisc(SAP.mod,term = "SorghumAccessionRed")

plot(SAP.can.g)
plot(SAP.can.t)


SAP.can.score <- SAP.can.t$scores
SAP.can.score.LN <- SAP.can.score[SAP.can.score$Treatment=='LN',]
SAP.can.score.HN <- SAP.can.score[SAP.can.score$Treatment=='HN',]
SAP.can.score.both <- merge(SAP.can.score.LN, SAP.can.score.HN, by="SorghumAccessionRed")

phe <- SAP.can.score.both[,c(1,3,5)]
colnames(phe) <- c("Taxa", "LN", "HN")
write.table(phe, "phe.txt", row.names = F, quote = F, sep = "\t")

library(rMVP)

MVP.Data(fileHMP = "SAP_imputed.hmp", filePhe = "phe.txt", sep.phe = "\t")
genotype <- attach.big.matrix("mvp.geno.desc")
phenotype <- read.table("mvp.phe", header = T)
map <- read.table("mvp.geno.map", header = T)
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    #nPC.MLM=3,
    nPC.FarmCPU = 4,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("FarmCPU")
  )
  gc()
}



















# omit omit omit 
# for canonical discriminant analysis
SAP.omit=na.omit(SAP)
geno.omit.count=as.data.frame(table(SAP.omit$SorghumAccessionRed))
geno.omit.keep = as.character(geno.omit.count[geno.omit.count$Freq==4,1])
SAP.o.k = SAP.omit[SAP.omit$SorghumAccessionRed %in% geno.omit.keep,]
SAP.o.k2=SAP.o.k[,c(3,8:26)]

pca.o <- princomp(SAP.o.k2[3:20],cor=T,scores = T)
#scores.o=cbind(SAP.avg[1:2],pca.o$scores)
plot(pca.o$scores,pch=SAP.o.k2$Treatment)

# CDA
SAP.o.k$Block = factor(as.numeric(SAP.o.k$Block<6.5))
SAP.o.k$Treatment = factor(SAP.o.k$Treatment)
SAP.o.k$SorghumAccessionRed = factor(SAP.o.k$SorghumAccessionRed)
SAP.o.cda = SAP.o.k[c(3,7,8:26)]
SAP.o.mod = lm(cbind(DaysToBloom,MedianLeafAngle,LeafAngleSDV,PaniclesPerPlot,
                     PanicleGrainWeight,EstimatedPlotYield,FlagLeafLength,
                     FlagLeafWidth,ExtantLeafNumber,PlantHeight,TillersPerPlant,
                     StemDiameterLower,StemDiameterUpper,RachisLength,
                     RachisDiameterLower,RachisDiameterUpper,PrimaryBranchNo,BranchInternodeLength)
               ~Treatment*SorghumAccessionRed+Block:Treatment,data=SAP.o.cda)
#SAP.o.mod = aov(as.formula(paste(paste(colnames(SAP.o.cda[-c(1:3)]),collapse = ","),"~Treatment*SorghumAccessionRed+Error(Treatment/Block)")),data=SAP.o.cda)
# SAP.o.mod = aov(cbind(DaysToBloom,MedianLeafAngle,LeafAngleSDV,PaniclesPerPlot,
#                       PanicleGrainWeight,EstimatedPlotYield,FlagLeafLength,
#                       FlagLeafWidth,ExtantLeafNumber,PlantHeight,TillersPerPlant,
#                       StemDiameterLower,StemDiameterUpper,RachisLength,
#                       RachisDiameterLower,RachisDiameterUpper,PrimaryBranchNo,BranchInternodeLength)
#                       ~Treatment*SorghumAccessionRed+Error(Block:Treatment),data=SAP.o.cda)
SAP.anova.out=Anova(SAP.o.mod,type='III',test="Wilks")
SAP.can.t=candisc(SAP.o.mod,term = "Treatment")
SAP.can.g=candisc(SAP.o.mod,term = "SorghumAccessionRed")
SAP.can.g1=candisc(SAP.o.mod,term = "SorghumAccessionRed",ndim=1)
SAP.can.ge=candisc(SAP.o.mod,term = "Treatment:SorghumAccessionRed")
SAP.can.ge1=candisc(SAP.o.mod,term = "Treatment:SorghumAccessionRed",ndim = 1)
SAP.can.t
SAP.can.t
SAP.can.g
plot(SAP.can.g)
plot(SAP.can.t)
plot(SAP.can.ge)
plot(SAP.can.ge1)
plot(SAP.can.g1)

library(MASS)
SAP.o.lda=SAP.o.cda[,-c(1,2)]
SAP.lda=lda(Treatment~.,data = SAP.o.lda)
SAP.lda



qqnorm(SAP.o.cda$DaysToBloom)
qqnorm(SAP.o.cda$MedianLeafAngle)
qqnorm(SAP.o.cda$LeafAngleSDV)
qqnorm(SAP.o.cda$PaniclesPerPlot)
qqnorm(SAP.o.cda$PanicleGrainWeight)
qqnorm(SAP.o.cda$EstimatedPlotYield)
qqnorm(SAP.o.cda$FlagLeafLength)
qqnorm(SAP.o.cda$FlagLeafWidth)
qqnorm(SAP.o.cda$ExtantLeafNumber)
qqnorm(SAP.o.cda$PlantHeight)
qqnorm(SAP.o.cda$TillersPerPlant)
qqnorm(SAP.o.cda$StemDiameterLower)
qqnorm(SAP.o.cda$StemDiameterUpper)
qqnorm(SAP.o.cda$RachisLength)
qqnorm(SAP.o.cda$RachisDiameterLower)
qqnorm(SAP.o.cda$RachisDiameterUpper)
qqnorm(SAP.o.cda$PrimaryBranchNo)
qqnorm(SAP.o.cda$BranchInternodeLength)

library(e1071)
trans.posskew=function(x){
  sk=skewness(x)
  if(sk<5 & sk >1){
    x=x^(1/3)
  }else if(sk>=5){
    x=log(x)
  }else{
    x
  }
}

SAP.ot.cda=SAP.o.cda
SAP.ot.cda[,-c(1:3)]=sapply(SAP.o.cda[,-c(1:3)], trans.posskew)
SAP.ot.mod = lm(cbind(DaysToBloom,MedianLeafAngle,LeafAngleSDV,PaniclesPerPlot,
                     PanicleGrainWeight,EstimatedPlotYield,FlagLeafLength,
                     FlagLeafWidth,ExtantLeafNumber,PlantHeight,TillersPerPlant,
                     StemDiameterLower,StemDiameterUpper,RachisLength,
                     RachisDiameterLower,RachisDiameterUpper,PrimaryBranchNo,BranchInternodeLength)
               ~Treatment*SorghumAccessionRed+Block:Treatment,data=SAP.ot.cda)
SAP.t.anova.out=Anova(SAP.ot.mod,type='III',test="Wilks")
SAP.t.can.t=candisc(SAP.ot.mod,term = "Treatment")
SAP.t.can.g=candisc(SAP.ot.mod,term = "SorghumAccessionRed")
SAP.t.can.ge=candisc(SAP.ot.mod,term = "Treatment:SorghumAccessionRed")
#SAP.t.can.t
#SAP.t.can.g
plot(SAP.t.can.g)
plot(SAP.t.can.t)
plot(SAP.t.can.ge)

