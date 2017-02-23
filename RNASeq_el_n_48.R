############################################################################
############################################################################
###                                                                      ###
### This file provides R code to reproduce the simulations               ###
### presented in Figures 3, 4 and 5 of                                   ###
###                                                                      ###
###                                                                      ###
###                                                                      ###
###                                                                      ###
############################################################################
############################################################################

############################################################################
### Install needed packages
###

#install.packages("doParallel")


library(DESeq)
library(DESeq2)
library(baySeq)
library(edgeR)
library(EBSeq)
library(NOISeq)
library(samr)
library(limma)
library(MASS)

## simulation set up

num<-1
total<-100
diffpro<-c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90)
ngene<-10000
ndiff <- ngene*diffpro
samplesize <- 48
trtsize <- samplesize/2
label <- c(rep(0, samplesize-trtsize), rep(1, trtsize))
#resample<-1000

## Generate simulation data
countsTableNS <- read.delim( "D:/Dongmei/Research/NGS/suppl_data/NeuralStemCellData.tab", row.names=1 )
#countsTableNS <- read.delim( "/scratch/dli28/suppl_data/NeuralStemCellData.tab", row.names=1 )
condsNS <- c( "GNS_dupe", "GNS", "GNS", "GNSL", "NS", "NS" )
NSmean <- apply(countsTableNS[, 5:6], 1, mean)
NScolsum <- apply(countsTableNS[, 5:6], 2, sum)
NSdep <- mean(NScolsum)
NSlambda <- NSmean/NSdep
NSvar <- apply(countsTableNS[, 5:6], 1, var)
theta <- NSmean^2/(NSvar-NSmean)
thetanew <- theta[is.finite(theta)]

equal <- TRUE

# Library size
if(equal){
	expected.lib.size <- rep(11e6, samplesize)
} else {
	expected.lib.size <- 20e6 * rep(c(1, 0.5), samplesize/2)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

set.seed(123)

### case 1: correlation coefficients among genes are 0

alpha<-0.05

totRejection_edgeR_exact<-FDR_edgeR_exact<-Power_edgeR_exact<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_edgeR_glm<-FDR_edgeR_glm<-Power_edgeR_glm<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_DESeq<-FDR_DESeq<-Power_DESeq<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_DESeq2<-FDR_DESeq2<-Power_DESeq2<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_baySeq<-FDR_baySeq<-Power_baySeq<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_EBSeq<-FDR_EBSeq<-Power_EBSeq<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_SAMSeq<-FDR_SAMSeq<-Power_SAMSeq<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_NOISeq<-FDR_NOISeq<-Power_NOISeq<-matrix(NA, nrow=length(diffpro), ncol=total)
totRejection_Voom<-FDR_Voom<-Power_Voom<-matrix(NA, nrow=length(diffpro), ncol=total)


for (j in 1:length(diffpro))
{
 num<-1
 fc <-c(rep(2, ndiff[j]/2), rep(1/2, ndiff[j]/2), rep(1, ngene-ndiff[j]))
  while (num<=total)
{
  control <- trt <- matrix(NA, nrow = ngene, ncol=trtsize)
  for (i in 1:samplesize/2)
  {
   mucontrol <- expected.lib.size[i]* NSlambda[1:ngene]
   mucontrol[mucontrol==0]<-0.5
   mutrt <- fc * mucontrol
   control[, i] <- rnegbin(ngene, mu = mucontrol, theta = thetanew[1:ngene])
   trt[, i] <- rnegbin(ngene, mu = mutrt, theta = thetanew[1:ngene])
  }
  group<-cbind(control, trt)
  group[is.na(group)] <- 0
  rownames(group) <- paste("Gene",1:ngene)
  colnames(group)<-c(paste("Control", 1:(samplesize-trtsize)), paste("Treatment", 1:trtsize))
  design <- cbind(Grp1=1,Grp2vs1=c(rep(0, (samplesize-trtsize)), rep(1, trtsize)))
  options(digit=3)

  ## Filter choice 1
  keep <- rowSums(group)>=10
  ndiffnew <- sum(rowSums(group[1:ndiff[j],])>=10)
  nkeep <- sum(keep)
  group2 <- group[keep,]

  ## Filter choice 2

  #group_fil <- filtered.data(group, factor=label, norm=FALSE, method = 3, cv.cutoff=100, cpm=1, p.adj = "fdr")
  # group_keep <- rownames(group_fil)
  # ## split string at non-digits
  #s <- strsplit(group_keep, "[^[:digit:]]")

  ## convert strings to numeric ("" become NA)
  #solution <- as.numeric(unlist(s))

  ## remove NA and duplicates
  #solution <- unique(solution[!is.na(solution)])
  # ndiffnew <- sum(solution <= ndiff[j])

  ## edgeR Exact Test 3.12.0
  
  edgeR.dgelist <- DGEList(counts = group2, group = factor(label))
  edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = "TMM")
  edgeR.disp <- estimateCommonDisp(edgeR.dgelist)
  edgeR.disp <- estimateTagwiseDisp(edgeR.disp)
  edgeR.test <- exactTest(edgeR.disp)
  edgeR.pvalues <- edgeR.test$table$PValue
  edgeR.adjp <- p.adjust(edgeR.pvalues, "BH")

  totRejection_edgeR_exact[j, num]<-sum(edgeR.adjp <= 0.05, na.rm=TRUE)
  FDR_edgeR_exact[j, num]<-sum(edgeR.adjp[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_edgeR_exact[j, num]
  Power_edgeR_exact[j, num]<-sum(edgeR.adjp[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  ## edgeR GLM Test 3.12.0

  edgeR.glm <- estimateGLMTrendedDisp(edgeR.dgelist, design)
  edgeR.glm <- estimateGLMTagwiseDisp(edgeR.glm, design)
  edgeR.glmf <- glmFit(edgeR.glm, design)
  edgeR.glmde <- glmLRT(edgeR.glmf)
  edgeR.glmpvalues <- edgeR.glmde$table$PValue
  edgeR.glmadjp <- p.adjust(edgeR.glmpvalues, "BH")

  totRejection_edgeR_glm[j, num]<-sum(edgeR.glmadjp <= 0.05, na.rm=TRUE)
  FDR_edgeR_glm[j, num]<-sum(edgeR.glmadjp[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_edgeR_glm[j, num]
  Power_edgeR_glm[j, num]<-sum(edgeR.glmadjp[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  ## DESeq  1.22.1

  group2.norm <- round(tmm(group2), 0)
  DESeq.cds <- newCountDataSet(countData = group2.norm, conditions = factor(label))
  DESeq.disp <- estimateSizeFactors(DESeq.cds)
  DESeq.disp <- estimateDispersions(DESeq.disp, sharingMode = "maximum", method = "pooled", fitType = "local")
  DESeq.test <- nbinomTest(DESeq.disp, "0", "1")
  DESeq.pval <- DESeq.test$pval
  DESeq.adjp <- p.adjust(DESeq.pval, "BH")

  totRejection_DESeq[j, num]<-sum(DESeq.adjp <= 0.05, na.rm=TRUE)
  FDR_DESeq[j, num]<-sum(DESeq.adjp[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_DESeq[j, num]
  Power_DESeq[j, num]<-sum(DESeq.adjp[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  ## DESeq2  1.10.1

  group.con <- data.frame(cbind(c(rep("control", samplesize/2), rep("Treatment", samplesize/2)), c(rep("single-read", samplesize))))
  colnames(group.con) <- c("condition", "type")
  rownames(group.con) <- c(paste("Control", 1:(samplesize-trtsize)), paste("Treatment", 1:trtsize))
  DESeq2.dds <- DESeqDataSetFromMatrix(countData = group2.norm, colData = group.con, design = ~condition)
  DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
  DESeq2.pval <- results(DESeq2.test)$pvalue
  DESeq2.adjp <- p.adjust(DESeq2.pval, "BH")

  totRejection_DESeq2[j, num]<-sum(DESeq2.adjp <= 0.05, na.rm=TRUE)
  FDR_DESeq2[j, num]<-sum(DESeq2.adjp[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_DESeq2[j, num]
  Power_DESeq2[j, num]<-sum(DESeq2.adjp[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  ## baySeq 2.4.1

  baySeq.CD <- new("countData", data = group2, replicates = label, groups =list(NDE = c(rep(1, samplesize)), DE = label))
  libsizes(baySeq.CD) <-getLibsizes(baySeq.CD, estimationType = "edgeR")
  baySeq.CD <- getPriors.NB(baySeq.CD, samplesize = 5000, equalDispersions = TRUE, estimation = "QL", cl = NULL)
  baySeq.CD <- getLikelihoods(baySeq.CD, cl = NULL, verbose=FALSE)
  baySeq.posteriors.DE <- exp(baySeq.CD@posteriors)[, 2]
  baySeq.table <- topCounts(baySeq.CD, group = "DE", FDR = 1)
  baySeq.FDR <- baySeq.table$FDR[match(rownames(group2), rownames(baySeq.table))]

  totRejection_baySeq[j, num]<-sum(baySeq.FDR <= 0.05, na.rm=TRUE)
  FDR_baySeq[j, num]<-sum(baySeq.FDR[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_baySeq[j, num]
  Power_baySeq[j, num]<-sum(baySeq.FDR[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  ## EBSeq  1.10.0

  sizes <- MedianNorm(group2.norm)
  quiet(EBSeq.out <- EBTest(Data = group2.norm, Conditions = factor(label), sizeFactors = sizes, maxround = 5))
  EBSeq.result <- GetDEResults(EBSeq.out, FDR = 0.05)
  EBSeq.status <- EBSeq.result$Status

  totRejection_EBSeq[j, num]<-sum(EBSeq.status == "DE", na.rm=TRUE)
  FDR_EBSeq[j, num]<-sum(EBSeq.status[(ndiffnew+1):dim(group2)[1]] == "DE", na.rm=TRUE)/totRejection_EBSeq[j, num]
  Power_EBSeq[j, num]<-sum(EBSeq.status[1:ndiffnew] == "DE", na.rm=TRUE)/ndiffnew

  ## SAMSeq   2.0
  label2<-c(rep(1, samplesize-trtsize), rep(2, trtsize))
  quiet(SAMSeq.test <- SAMseq(group2.norm, label2, resp.type = "Two class unpaired", geneid = rownames(group2.norm), genenames = rownames(group2.norm), nperms = 100, nresamp = 20, fdr.output = 0.05))
  SAMSeq.result <- rbind(SAMSeq.test$siggenes.table$genes.up, SAMSeq.test$siggenes.table$genes.lo)
  SAMSeq.s <- strsplit(SAMSeq.result[, 1], "[^[:digit:]]")
  SAMSeq.solution <- as.numeric(unlist(SAMSeq.s))
  SAMSeq.solution <- unique(SAMSeq.solution[!is.na(SAMSeq.solution)])

  totRejection_SAMSeq[j, num]<-length(SAMSeq.solution)
  FDR_SAMSeq[j, num]<-sum(SAMSeq.solution > ndiff[j])/totRejection_SAMSeq[j, num]
  Power_SAMSeq[j, num]<-sum(SAMSeq.solution <= ndiff[j])/ndiffnew

  ## NOISeq    2.14.1
  NOISeq.data <- readData(group2, factors = group.con)
  quiet(NOISeq.result <- noiseqbio(NOISeq.data, norm = "tmm", factor = "condition", lc = 1, filter = 0))
  NOISeq.DE <- degenes(NOISeq.result, q = 0.95, M = NULL)
  NOISeq.s <- strsplit(rownames(NOISeq.DE), "[^[:digit:]]")
  NOISeq.solution <- as.numeric(unlist(NOISeq.s))
  NOISeq.solution <- unique(NOISeq.solution[!is.na(NOISeq.solution)])

  totRejection_NOISeq[j, num]<-length(NOISeq.solution)
  FDR_NOISeq[j, num]<-sum(NOISeq.solution > ndiff[j])/totRejection_NOISeq[j, num]
  Power_NOISeq[j, num]<-sum(NOISeq.solution <= ndiff[j])/ndiffnew

  ##voom + limma  3.26.8

  Voom.data <- voom(group2.norm, design = design)
  Voom.limma <- lmFit(Voom.data, design, na.rm=TRUE)
  Voom.ebayes<- eBayes(Voom.limma)
  Voom.pval<-Voom.ebayes$p.value[, 2]
  Voom.adjp<-p.adjust(Voom.pval, method = "BH")

  totRejection_Voom[j, num]<-sum(Voom.adjp <= 0.05, na.rm=TRUE)
  FDR_Voom[j, num]<-sum(Voom.adjp[(ndiffnew+1):dim(group2)[1]] <= 0.05, na.rm=TRUE)/totRejection_Voom[j, num]
  Power_Voom[j, num]<-sum(Voom.adjp[1:ndiffnew] <= 0.05, na.rm=TRUE)/ndiffnew

  num<-num+1
 }
}

meantotRejection_edgeR_exact <- vartotRejection_edgeR_exact <- meanFDR_edgeR_exact <- meanPower_edgeR_exact <- rep(NA, length(diffpro))
meantotRejection_edgeR_glm <- vartotRejection_edgeR_glm <- meanFDR_edgeR_glm <- meanPower_edgeR_glm <- rep(NA, length(diffpro))
meantotRejection_DESeq <- vartotRejection_DESeq <- meanFDR_DESeq <- meanPower_DESeq <- rep(NA, length(diffpro))
meantotRejection_DESeq2 <- vartotRejection_DESeq2 <- meanFDR_DESeq2 <- meanPower_DESeq2 <- rep(NA, length(diffpro))
meantotRejection_baySeq <- vartotRejection_baySeq <- meanFDR_baySeq <- meanPower_baySeq <- rep(NA, length(diffpro))
meantotRejection_EBSeq <- vartotRejection_EBSeq <- meanFDR_EBSeq <- meanPower_EBSeq <- rep(NA, length(diffpro))
meantotRejection_SAMSeq <- vartotRejection_SAMSeq <- meanFDR_SAMSeq <- meanPower_SAMSeq <- rep(NA, length(diffpro))
meantotRejection_NOISeq <- vartotRejection_NOISeq <- meanFDR_NOISeq <- meanPower_NOISeq <- rep(NA, length(diffpro))
meantotRejection_Voom <- vartotRejection_Voom <- meanFDR_Voom <- meanPower_Voom <- rep(NA, length(diffpro))

for (i in 1:length(diffpro))
{
  meanFDR_edgeR_exact[i]<-mean(FDR_edgeR_exact[i, ], na.rm=TRUE)
  meanPower_edgeR_exact[i]<-mean(Power_edgeR_exact[i, ], na.rm=TRUE)
  meantotRejection_edgeR_exact[i]<-mean(totRejection_edgeR_exact[i, ], na.rm=TRUE)
  vartotRejection_edgeR_exact[i]<-sqrt(var(totRejection_edgeR_exact[i, ], na.rm=TRUE))

  meanFDR_edgeR_glm[i]<-mean(FDR_edgeR_glm[i, ], na.rm=TRUE)
  meanPower_edgeR_glm[i]<-mean(Power_edgeR_glm[i, ], na.rm=TRUE)
  meantotRejection_edgeR_glm[i]<-mean(totRejection_edgeR_glm[i, ], na.rm=TRUE)
  vartotRejection_edgeR_glm[i]<-sqrt(var(totRejection_edgeR_glm[i, ], na.rm=TRUE))

  meanFDR_DESeq[i]<-mean(FDR_DESeq[i, ], na.rm=TRUE)
  meanPower_DESeq[i]<-mean(Power_DESeq[i, ], na.rm=TRUE)
  meantotRejection_DESeq[i]<-mean(totRejection_DESeq[i, ], na.rm=TRUE)
  vartotRejection_DESeq[i]<-sqrt(var(totRejection_DESeq[i, ], na.rm=TRUE))

  meanFDR_DESeq2[i]<-mean(FDR_DESeq2[i, ], na.rm=TRUE)
  meanPower_DESeq2[i]<-mean(Power_DESeq2[i, ], na.rm=TRUE)
  meantotRejection_DESeq2[i]<-mean(totRejection_DESeq2[i, ], na.rm=TRUE)
  vartotRejection_DESeq2[i]<-sqrt(var(totRejection_DESeq2[i, ], na.rm=TRUE))

  meanFDR_baySeq[i]<-mean(FDR_baySeq[i, ], na.rm=TRUE)
  meanPower_baySeq[i]<-mean(Power_baySeq[i, ], na.rm=TRUE)
  meantotRejection_baySeq[i]<-mean(totRejection_baySeq[i, ], na.rm=TRUE)
  vartotRejection_baySeq[i]<-sqrt(var(totRejection_baySeq[i, ], na.rm=TRUE))

  meanFDR_EBSeq[i]<-mean(FDR_EBSeq[i, ], na.rm=TRUE)
  meanPower_EBSeq[i]<-mean(Power_EBSeq[i, ], na.rm=TRUE)
  meantotRejection_EBSeq[i]<-mean(totRejection_EBSeq[i, ], na.rm=TRUE)
  vartotRejection_EBSeq[i]<-sqrt(var(totRejection_EBSeq[i, ], na.rm=TRUE))

  meanFDR_SAMSeq[i]<-mean(FDR_SAMSeq[i, ], na.rm=TRUE)
  meanPower_SAMSeq[i]<-mean(Power_SAMSeq[i, ], na.rm=TRUE)
  meantotRejection_SAMSeq[i]<-mean(totRejection_SAMSeq[i, ], na.rm=TRUE)
  vartotRejection_SAMSeq[i]<-sqrt(var(totRejection_SAMSeq[i, ], na.rm=TRUE))

  meanFDR_NOISeq[i]<-mean(FDR_NOISeq[i, ], na.rm=TRUE)
  meanPower_NOISeq[i]<-mean(Power_NOISeq[i, ], na.rm=TRUE)
  meantotRejection_NOISeq[i]<-mean(totRejection_NOISeq[i, ], na.rm=TRUE)
  vartotRejection_NOISeq[i]<-sqrt(var(totRejection_NOISeq[i, ], na.rm=TRUE))

  meanFDR_Voom[i]<-mean(FDR_Voom[i, ], na.rm=TRUE)
  meanPower_Voom[i]<-mean(Power_Voom[i, ], na.rm=TRUE)
  meantotRejection_Voom[i]<-mean(totRejection_Voom[i, ], na.rm=TRUE)
  vartotRejection_Voom[i]<-sqrt(var(totRejection_Voom[i, ], na.rm=TRUE))
}

print(cbind(meanFDR_edgeR_exact, meanFDR_edgeR_glm, meanFDR_DESeq, meanFDR_DESeq2,  meanFDR_baySeq, meanFDR_EBSeq, meanFDR_SAMSeq, meanFDR_NOISeq, meanFDR_Voom))
print(cbind(meanPower_edgeR_exact, meanPower_edgeR_glm, meanPower_DESeq, meanPower_DESeq2,  meanPower_baySeq, meanPower_EBSeq, meanPower_SAMSeq, meanPower_NOISeq, meanPower_Voom))
print(cbind(meantotRejection_edgeR_exact, meantotRejection_edgeR_glm, meantotRejection_DESeq, meantotRejection_DESeq2,  meantotRejection_baySeq, meantotRejection_EBSeq, meantotRejection_SAMSeq, meantotRejection_NOISeq, meantotRejection_Voom ))
print(cbind(vartotRejection_edgeR_exact, vartotRejection_edgeR_glm, vartotRejection_DESeq, vartotRejection_DESeq2,  vartotRejection_baySeq, vartotRejection_EBSeq, vartotRejection_SAMSeq, vartotRejection_NOISeq, vartotRejection_Voom))

meanFDR<-cbind(meanFDR_edgeR_exact, meanFDR_edgeR_glm, meanFDR_DESeq, meanFDR_DESeq2,  meanFDR_baySeq, meanFDR_EBSeq, meanFDR_SAMSeq, meanFDR_NOISeq, meanFDR_Voom)
meanPower<-cbind(meanPower_edgeR_exact, meanPower_edgeR_glm, meanPower_DESeq, meanPower_DESeq2,  meanPower_baySeq, meanPower_EBSeq, meanPower_SAMSeq, meanPower_NOISeq, meanPower_Voom)
meantotreject<-cbind(meantotRejection_edgeR_exact, meantotRejection_edgeR_glm, meantotRejection_DESeq, meantotRejection_DESeq2,  meantotRejection_baySeq, meantotRejection_EBSeq, meantotRejection_SAMSeq, meantotRejection_NOISeq, meantotRejection_Voom )
stdreject<-cbind(vartotRejection_edgeR_exact, vartotRejection_edgeR_glm, vartotRejection_DESeq, vartotRejection_DESeq2,  vartotRejection_baySeq, vartotRejection_EBSeq, vartotRejection_SAMSeq, vartotRejection_NOISeq, vartotRejection_Voom)

#write.csv(meanFDR, file="meanFDR_el_n_6.csv")
#write.csv(meanPower, file="meanPower_el_n_6.csv")
#write.csv(meantotreject, file="meantotreject_el_n_6.csv")
#write.csv(stdreject, file="stdreject_el_n_6.csv")

write.csv(meanFDR, file="D:/Dongmei/Research/NGS/RNASeq_simulation/meanFDR_el_n_24.csv")
write.csv(meanPower, file="D:/Dongmei/Research/NGS/RNASeq_simulation/meanPower_el_n_24.csv")
write.csv(meantotreject, file="D:/Dongmei/Research/NGS/RNASeq_simulation/meantotreject_el_n_24.csv")
write.csv(stdreject, file="D:/Dongmei/Research/NGS/RNASeq_simulation/stdreject_el_n_24.csv")


colnames(meanFDR) <- colnames(meanPower) <- colnames(meantotreject) <- colnames(stdreject) <- c("edgeR Exact", "edgeR GLM", "DESeq", "DESeq2", "baySeq", "EBSeq", "SAMSeq", "NOISeq", "Voom")

pdf(file="RNA_Seq_comparison_el_n_24.pdf", width=4, height=4)
#par(mfrow=c(2,2), cex=0.7, mar=c(5, 5, 2.5, 2.5))
par(mfrow=c(2, 2), mar=c(8,6,2,2)+0.1, mgp=c(4, 1, 0), cex=0.7)
boxplot(meanFDR, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="FDR", las=2)
boxplot(meanPower, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="Power", las=2)
boxplot(meantotreject, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="Total discoveries", las=2)
boxplot(stdreject, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="SD of total discoveries", las=2)
dev.off()

postscript(file="RNA_Seq_comparison_el_n_24.eps", width=4, height=4)
par(mfrow=c(2,2), cex=0.7, mar=c(4.5, 4.5, 1.5, 0.7))
boxplot(meanFDR, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="FDR", las=2)
boxplot(meanPower, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="Power", las=2)
boxplot(meantotreject, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="Total discoveries", las=2)
boxplot(stdreject, col=c("blue","red", "darkgreen", "cyan", "grey", "orange", "yellowgreen", "purple", "brown4"), xlab=NULL, ylab="SD of total discoveries", las=2)
dev.off()


pdf(file="RNA_Seq_comparison_el_n_24_lineplot.pdf", width=4, height=4)
par(mfrow=c(2,2), cex=0.7, mgp=c(3, 1, 0), mar=c(4.5, 4.5, 1.5, 0.7))
g_range<-range(meanFDR, na.rm=TRUE)
plot(diffpro, meanFDR_edgeR_exact, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, meanFDR_edgeR_glm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, meanFDR_DESeq, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, meanFDR_DESeq2, type="l", lty=4, lwd=2, col="black")
lines(diffpro, meanFDR_baySeq, type="l", lty=5, lwd=2, col="grey")
lines(diffpro, meanFDR_EBSeq, type="l", lty=6, lwd=2, col="orange")
lines(diffpro, meanFDR_SAMSeq, type="l", lty=1, lwd=2, col="yellowgreen")
lines(diffpro, meanFDR_NOISeq, type="l", lty=2, lwd=2, col="purple")
lines(diffpro, meanFDR_Voom, type="l", lty=3, lwd=2, col="brown4")
title(xlab=expression(pi[1]), ylab="FDR")
legend(0.6, g_range[2], c("EdgeR Exact","EdgeR GLM", "DESeq", "DESeq2", "baySeq", "EBSeq", "SAMSeq", "NOISeq", "Voom"), cex=0.8,
col=c("blue","red", "darkgreen", "black", "grey", "orange", "yellowgreen", "purple", "brown4"), lwd=2, lty=c(1:6, 1:3), bty="n")
#text(1, g_range[2], "A")


g_range<-range(meanPower, na.rm=TRUE)
plot(diffpro, meanPower_edgeR_exact, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, meanPower_edgeR_glm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, meanPower_DESeq, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, meanPower_DESeq2, type="l", lty=4, lwd=2, col="black")
lines(diffpro, meanPower_baySeq, type="l", lty=5, lwd=2, col="grey")
lines(diffpro, meanPower_EBSeq, type="l", lty=6, lwd=2, col="orange")
lines(diffpro, meanPower_SAMSeq, type="l", lty=1, lwd=2, col="yellowgreen")
lines(diffpro, meanPower_NOISeq, type="l", lty=2, lwd=2, col="purple")
lines(diffpro, meanPower_Voom, type="l", lty=3, lwd=2, col="brown4")
title(xlab=expression(pi[1]))
title(ylab="Power")
#legend(0.01, g_range[2], c("EdgeR Exact","EdgeR GLM", "DESeq", "DESeq2", "baySeq", "EBSeq", "SAMSeq", "NOISeq", "Voom"), cex=0.8,
#col=c("blue","red", "darkgreen", "black", "grey", "orange", "yellowgreen", "purple", "brown4"), lwd=2, lty=c(1:6, 1:3), bty="n")
#text(1, g_range[2], "B")


g_range<-range(meantotreject, na.rm=TRUE)
plot(diffpro, meantotRejection_edgeR_exact, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, meantotRejection_edgeR_glm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, meantotRejection_DESeq, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, meantotRejection_DESeq2, type="l", lty=4, lwd=2, col="black")
lines(diffpro, meantotRejection_baySeq, type="l", lty=5, lwd=2, col="grey")
lines(diffpro, meantotRejection_EBSeq, type="l", lty=6, lwd=2, col="orange")
lines(diffpro, meantotRejection_SAMSeq, type="l", lty=1, lwd=2, col="yellowgreen")
lines(diffpro, meantotRejection_NOISeq, type="l", lty=2, lwd=2, col="purple")
lines(diffpro, meantotRejection_Voom, type="l", lty=3, lwd=2, col="brown4")
title(xlab=expression(pi[1]))
title(ylab="Mean of total discoveries")
#legend(0.01, g_range[2], c("EdgeR Exact","EdgeR GLM", "DESeq", "DESeq2", "baySeq", "EBSeq", "SAMSeq", "NOISeq", "Voom"), cex=0.8,
#col=c("blue","red", "darkgreen", "black", "grey", "orange", "yellowgreen", "purple", "brown4"), lwd=2, lty=c(1:6, 1:3), bty="n")
#text(1, g_range[2], "C")

g_range<-range(stdreject, na.rm=TRUE)
plot(diffpro, vartotRejection_edgeR_exact, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, vartotRejection_edgeR_glm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, vartotRejection_DESeq, type="l", lty=3, lwd=2, col="dark green")
lines(diffpro, vartotRejection_DESeq2, type="l", lty=4, lwd=2, col="black")
lines(diffpro, vartotRejection_baySeq, type="l", lty=5, lwd=2, col="grey")
lines(diffpro, vartotRejection_EBSeq, type="l", lty=6, lwd=2, col="orange")
lines(diffpro, vartotRejection_SAMSeq, type="l", lty=1, lwd=2, col="yellowgreen")
lines(diffpro, vartotRejection_NOISeq, type="l", lty=2, lwd=2, col="purple")
lines(diffpro, vartotRejection_Voom, type="l", lty=3, lwd=2, col="brown4")
title(xlab=expression(pi[1]))
title(ylab="SD of total discoveries")
#legend(0.6, g_range[2], c("EdgeR Exact","EdgeR GLM", "DESeq", "DESeq2", "baySeq", "EBSeq", "SAMSeq", "NOISeq", "Voom"), cex=0.8,
#col=c("blue","red", "darkgreen", "black", "grey", "orange", "yellowgreen", "purple", "brown4"), lwd=2, lty=c(1:6, 1:3), bty="n")
#text(1, g_range[2], "D")

dev.off()

