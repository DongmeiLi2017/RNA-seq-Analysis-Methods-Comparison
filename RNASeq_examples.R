source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite("baySeq")
biocLite("EBSeq")
biocLite("NOISeq")
biocLite("samr")

set.seed(123)

## Import required library
library(DESeq)
library(DESeq2)
library(baySeq)
library(edgeR)
library(EBSeq)
library(NOISeq)
library(samr)
library(limma)
library(MASS)
library(VennDiagram)

#####################################################################################################################
#
# RNASeq comparison paper real data examples  Bottomly data
#
#####################################################################################################################

## Read the RNA-Seq data into R

bottomly_count <- read.table("D:/Dongmei/Research/Multiple_testing_procedures/new_MTP/new_MTP_examples/bottomly_count_table.txt", header=TRUE)
bottomly_id <- read.table("D:/Dongmei/Research/Multiple_testing_procedures/new_MTP/new_MTP_examples/bottomly_id_table.txt", header=TRUE)

#> dim(bottomly_count)
#[1] 36536    22

keep <- rowSums(bottomly_count[, -1])>=10
nkeep <- sum(keep)
#> nkeep
#[1] 11870
bottomly_countnew <- bottomly_count[keep,]
label <- c(rep(0, 10), rep(1, 11))
design <- cbind(Grp1=1,Grp2vs1=label)

## edgeR Exact Test 3.12.0

edgeR.dgelist <- DGEList(counts = bottomly_countnew[, -1], group = factor(label))
edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.disp <- estimateCommonDisp(edgeR.dgelist)
edgeR.disp <- estimateTagwiseDisp(edgeR.disp)
edgeR.test <- exactTest(edgeR.disp)
edgeR.pvalues <- edgeR.test$table$PValue
edgeR.adjp <- p.adjust(edgeR.pvalues, "BH")
sum(edgeR.adjp <=0.05)
#[1] 1235
edgeR.bottomly <- cbind(bottomly_countnew, edgeR.test$table, edgeR.adjp)
edgeR.bottomly.select <- edgeR.bottomly[which(edgeR.adjp <=0.05), ]

edgeR.sig.gene <- as.numeric(rownames(bottomly.norm)[which(edgeR.adjp <=0.05)])

## edgeR GLM Test 3.12.0

edgeR.glm <- estimateGLMTrendedDisp(edgeR.dgelist, design)
edgeR.glm <- estimateGLMTagwiseDisp(edgeR.glm, design)
edgeR.glmf <- glmFit(edgeR.glm, design)
edgeR.glmde <- glmLRT(edgeR.glmf)
edgeR.glmpvalues <- edgeR.glmde$table$PValue
edgeR.glmadjp <- p.adjust(edgeR.glmpvalues, "BH")
sum(edgeR.glmadjp <=0.05)
#[1] 1262
edgeR.glm.bottomly <- cbind(bottomly_countnew, edgeR.glmde$table, edgeR.glmadjp)
edgeR.glm.bottomly.select <- edgeR.glm.bottomly[which(edgeR.glmadjp <=0.05), ]

edgeR.glm.sig.gene <- as.numeric(rownames(bottomly.norm)[which(edgeR.glmadjp <=0.05)])

## DESeq  1.22.1

bottomly.norm <- round(edgeR.dgelist$counts, 0)
DESeq.cds <- newCountDataSet(countData = bottomly.norm, conditions = factor(label))
DESeq.disp <- estimateSizeFactors(DESeq.cds)
DESeq.disp <- estimateDispersions(DESeq.disp, sharingMode = "maximum", method = "pooled", fitType = "local")
DESeq.test <- nbinomTest(DESeq.disp, "0", "1")
DESeq.pval <- DESeq.test$pval
DESeq.adjp <- p.adjust(DESeq.pval, "BH")
sum(DESeq.adjp <=0.05)
#[1] 598
DESeq.bottomly <- cbind(bottomly_countnew, DESeq.test, DESeq.adjp)
DESeq.bottomly.select <- DESeq.bottomly[which(DESeq.adjp <=0.05), ]

DESeq.sig.gene <- as.numeric(rownames(bottomly.norm)[which(DESeq.adjp <= 0.05)])

## DESeq2  1.10.1

group.con <- data.frame(cbind(c(rep("C57BL", 10), rep("DBA", 11)), c(rep("single-read", 21))))
colnames(group.con) <- c("condition", "type")
rownames(group.con) <- c(paste("C57BL", 1:10), paste("DBA", 1:11))
colnames(bottomly.norm) <- rownames(group.con)
DESeq2.dds <- DESeqDataSetFromMatrix(countData = bottomly.norm, colData = group.con, design = ~condition)
DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
DESeq2.pval <- results(DESeq2.test)$pvalue
DESeq2.adjp <- p.adjust(DESeq2.pval, "BH")
sum(DESeq2.adjp <= 0.05, na.rm=TRUE)
#[1] 1314  1325
DESeq2.bottomly <- cbind(bottomly_countnew, results(DESeq2.test), DESeq2.adjp)
DESeq2.bottomly.select <- DESeq2.bottomly[which(DESeq2.adjp <= 0.05), ]

DESeq2.sig.gene <- as.numeric(rownames(bottomly.norm)[which(DESeq2.adjp <= 0.05)])

## baySeq 2.4.1

baySeq.CD <- new("countData", data = bottomly.norm, replicates = label, groups =list(NDE = c(rep(1, 21)), DE = label))
libsizes(baySeq.CD) <-getLibsizes(baySeq.CD, estimationType = "edgeR")
baySeq.CD <- getPriors.NB(baySeq.CD, samplesize = 5000, equalDispersions = TRUE, estimation = "QL", cl = NULL)
baySeq.CD <- getLikelihoods(baySeq.CD, cl = NULL, verbose=FALSE)
baySeq.posteriors.DE <- exp(baySeq.CD@posteriors)[, 2]
baySeq.table <- topCounts(baySeq.CD, group = "DE", FDR = 1)
baySeq.FDR <- baySeq.table$FDR[match(rownames(bottomly.norm), rownames(baySeq.table))]
sum(baySeq.FDR <= 0.05, na.rm=TRUE)
# [1] 521 523 519
baySeq.bottomly <- cbind(bottomly_countnew, baySeq.table)
baySeq.bottomly.select <- baySeq.bottomly[which(baySeq.FDR <= 0.05), ]

baySeq.sig.gene <- as.numeric(rownames(bottomly.norm)[which(baySeq.FDR <= 0.05)])

## EBSeq  1.10.0

sizes <- MedianNorm(bottomly.norm)
EBSeq.out <- EBTest(Data = bottomly.norm, Conditions = factor(label), sizeFactors = sizes, maxround = 5)
EBSeq.result <- GetDEResults(EBSeq.out, FDR = 0.05)
EBSeq.status <- EBSeq.result$Status
sum(EBSeq.status == "DE", na.rm=TRUE)
#[1] 378
EBSeq.adjp <- EBSeq.result$PPMat[, 1]
sum(EBSeq.adjp <= 0.05, na.rm=TRUE)
#[1] 378
EBSeq.bottomly <- cbind(bottomly_countnew, EBSeq.out$PPDE, EBSeq.out$f0, EBSeq.out$f1, EBSeq.status)
EBSeq.bottomly.select <- EBSeq.bottomly[which(EBSeq.status == "DE"), ]

EBSeq.sig.gene <- as.numeric(rownames(bottomly.norm)[which(EBSeq.status == "DE")])

## SAMSeq   2.0
label2<-c(rep(1, 10), rep(2, 11))
SAMSeq.test <- SAMseq(bottomly.norm, label2, resp.type = "Two class unpaired", geneid = rownames(bottomly.norm), genenames = rownames(bottomly.norm), nperms = 100, nresamp = 20, fdr.output = 0.05)
SAMSeq.result <- rbind(SAMSeq.test$siggenes.table$genes.up, SAMSeq.test$siggenes.table$genes.lo)
SAMSeq.s <- strsplit(SAMSeq.result[, 1], "[^[:digit:]]")
SAMSeq.solution <- as.numeric(unlist(SAMSeq.s))
SAMSeq.solution <- unique(SAMSeq.solution[!is.na(SAMSeq.solution)])
length(SAMSeq.solution)
#[1] 1968    1775  1933  1727   1779
SAMSeq.bottomly.select <- bottomly_countnew[SAMSeq.solution, ]

SAMSeq.sig.gene <- SAMSeq.solution

SAMSeq.adjp <- rep(1, length(EBSeq.adjp))
SAMSeq.adjp[SAMSeq.solution] <- 0

## NOISeq    2.14.1
NOISeq.data <- readData(bottomly.norm, factors = group.con)
NOISeq.result <- noiseqbio(NOISeq.data, norm = "tmm", factor = "condition", lc = 1, filter = 0)
NOISeq.DE <- degenes(NOISeq.result, q = 0.95, M = NULL)
NOISeq.s <- strsplit(rownames(NOISeq.DE), "[^[:digit:]]")
NOISeq.solution <- as.numeric(unlist(NOISeq.s))
NOISeq.solution <- unique(NOISeq.solution[!is.na(NOISeq.solution)])
length(NOISeq.solution)
#[1] 409

NOISeq.adjp <- rep(1, length(EBSeq.adjp))
NOISeq.adjp[NOISeq.solution] <- 0

#NOISeq.bottomly <- cbind(bottomly_countnew, NOISeq.result)
#NOISeq.bottomly.select <- EBSeq.bottomly[which(EBSeq.status == "DE"), ]

NOISeq.bottomly.select <- bottomly_countnew[NOISeq.solution, ]

NOISeq.sig.gene <- NOISeq.solution

##voom + limma  3.26.8

Voom.data <- voom(bottomly.norm, design = design)
Voom.limma <- lmFit(Voom.data, design, na.rm=TRUE)
Voom.ebayes<- eBayes(Voom.limma)
Voom.pval<-Voom.ebayes$p.value[, 2]
Voom.adjp<-p.adjust(Voom.pval, method = "BH")
sum(Voom.adjp <= 0.05, na.rm=TRUE)
#[1] 1037

Voom.bottomly <- cbind(bottomly_countnew, Voom.ebayes, Voom.adjp)
Voom.bottomly.select <- Voom.bottomly[which(Voom.adjp <= 0.05), ]

Voom.sig.gene <- as.numeric(rownames(bottomly.norm)[which(Voom.adjp <= 0.05)])


## Check correlation of adjusted p-values from different methods

adjp_all <- cbind(edgeR.adjp, edgeR.glmadjp, DESeq.adjp, DESeq2.adjp, baySeq.FDR, EBSeq.adjp, Voom.adjp, SAMSeq.adjp)
cor(adjp_all, use = "pairwise.complete.obs", method = "pearson")
#              edgeR.adjp edgeR.glmadjp DESeq.adjp DESeq2.adjp baySeq.FDR EBSeq.adjp Voom.adjp SAMSeq.adjp
#edgeR.adjp     1.0000000     0.9925868  0.8947244   0.9886325  0.7554290  0.5662760 0.9442558   0.6420589
#edgeR.glmadjp  0.9925868     1.0000000  0.8987626   0.9875338  0.7578535  0.5621737 0.9493447   0.6320304
#DESeq.adjp     0.8947244     0.8987626  1.0000000   0.8919513  0.9046539  0.7487988 0.8659809   0.7535920
#DESeq2.adjp    0.9886325     0.9875338  0.8919513   1.0000000  0.7475949  0.5702753 0.9394631   0.6556711
#baySeq.FDR     0.7554290     0.7578535  0.9046539   0.7475949  1.0000000  0.8297092 0.7344083   0.7649590
#EBSeq.adjp     0.5662760     0.5621737  0.7487988   0.5702753  0.8297092  1.0000000 0.5515724   0.6998865
#Voom.adjp      0.9442558     0.9493447  0.8659809   0.9394631  0.7344083  0.5515724 1.0000000   0.6261583
#SAMSeq.adjp    0.6420589     0.6320304  0.7535920   0.6556711  0.7649590  0.6998865 0.6261583   1.0000000

cor(adjp_all, use = "pairwise.complete.obs", method = "kendall")
#              edgeR.adjp edgeR.glmadjp DESeq.adjp DESeq2.adjp baySeq.FDR EBSeq.adjp Voom.adjp SAMSeq.adjp
#edgeR.adjp     1.0000000     0.9321778  0.7952570   0.9116007  0.5327397  0.4033518 0.7800348   0.4847668
#edgeR.glmadjp  0.9321778     1.0000000  0.8262929   0.9158791  0.5619612  0.4064123 0.7901427   0.4810586
#DESeq.adjp     0.7952570     0.8262929  1.0000000   0.8031200  0.5634135  0.4447182 0.7400546   0.4738940
#DESeq2.adjp    0.9116007     0.9158791  0.8031200   1.0000000  0.5256357  0.4085858 0.7642554   0.4910190
#baySeq.FDR     0.5327397     0.5619612  0.5634135   0.5256357  1.0000000  0.5335100 0.5232537   0.4413227
#EBSeq.adjp     0.4033518     0.4064123  0.4447182   0.4085858  0.5335100  1.0000000 0.3719779   0.3620765
#Voom.adjp      0.7800348     0.7901427  0.7400546   0.7642554  0.5232537  0.3719779 1.0000000   0.4807136
#SAMSeq.adjp    0.4847668     0.4810586  0.4738940   0.4910190  0.4413227  0.3620765 0.4807136   1.0000000

cor(adjp_all, use = "pairwise.complete.obs", method = "spearman")
#              edgeR.adjp edgeR.glmadjp DESeq.adjp DESeq2.adjp baySeq.FDR EBSeq.adjp Voom.adjp SAMSeq.adjp
#edgeR.adjp     1.0000000     0.9904175  0.9392125   0.9875060  0.7044914  0.5469865 0.9305726   0.5922510
#edgeR.glmadjp  0.9904175     1.0000000  0.9555208   0.9883360  0.7343896  0.5472660 0.9378934   0.5890615
#DESeq.adjp     0.9392125     0.9555208  1.0000000   0.9454573  0.7329498  0.5865164 0.9054984   0.5735923
#DESeq2.adjp    0.9875060     0.9883360  0.9454573   1.0000000  0.6995446  0.5563886 0.9245374   0.6012646
#baySeq.FDR     0.7044914     0.7343896  0.7329498   0.6995446  1.0000000  0.6630922 0.6960723   0.5404850
#EBSeq.adjp     0.5469865     0.5472660  0.5865164   0.5563886  0.6630922  1.0000000 0.5081914   0.4434207
#Voom.adjp      0.9305726     0.9378934  0.9054984   0.9245374  0.6960723  0.5081914 1.0000000   0.5886407
#SAMSeq.adjp    0.5922510     0.5890615  0.5735923   0.6012646  0.5404850  0.4434207 0.5886407   1.0000000


write.csv(edgeR.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_edgeR_select.csv")
write.csv(edgeR.glm.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_edgeR_glm_select.csv")
write.csv(DESeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_DESeq_select.csv")
write.csv(DESeq2.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_DESeq2_select.csv")
write.csv(baySeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_baySeq_select.csv")
write.csv(EBSeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_EBSeq_select.csv")
write.csv(SAMSeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_SAMSeq_select.csv")
write.csv(NOISeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_NOISeq_select.csv")
write.csv(Voom.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_Voom_select.csv")


select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "edgeR_glm" = which(edgeR.glmadjp <=0.05), "DESeq" = which(DESeq.adjp <= 0.05),
                                   "DESeq2" = which(DESeq2.adjp <= 0.05), "baySeq" = which(baySeq.FDR <= 0.05), "EBSeq"= which(EBSeq.status == "DE"),
                                   "SAMSeq" = SAMSeq.solution, "NOISeq" = NOISeq.solution, "Voom" = which(Voom.adjp <= 0.05)))

select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "edgeR_glm" = which(edgeR.glmadjp <=0.05), "DESeq" = which(DESeq.adjp <= 0.05),
                                   "DESeq2" = which(DESeq2.adjp <= 0.05), "baySeq" = which(baySeq.FDR <= 0.05), "EBSeq"= which(EBSeq.status == "DE"),
                                   "SAMSeq" = SAMSeq.solution,  "Voom" = which(Voom.adjp <= 0.05)))

select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "Voom" = which(Voom.adjp <=0.05), "DESeq2" = which(DESeq2.adjp <= 0.05),
                                          "EBSeq"= which(EBSeq.status == "DE"), "SAMSeq" = SAMSeq.solution))

## venn diagram plot
                    
venn.plot <- venn.diagram(x=list("edgeR Exact" = edgeR.sig.gene, "Voom" = Voom.sig.gene, "DESeq2" = DESeq2.sig.gene,
                    "EBSeq"= EBSeq.sig.gene, "SAMSeq" = SAMSeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    alpha = 0.50,
                    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                            1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)
                    
                    
venn.plot <- venn.diagram(x=list("edgeR Exact" = edgeR.sig.gene, "Voom" = Voom.sig.gene, "DESeq2" = DESeq2.sig.gene,
                    "EBSeq"= EBSeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot_without_SAMSeq.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    alpha = 0.50,
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)
                    

venn.plot <- venn.diagram(x=list("DESeq2" = DESeq2.sig.gene, "DESeq" = DESeq.sig.gene, "baySeq" = baySeq.sig.gene, 
                    "SAMSeq" = SAMSeq.sig.gene, "NOISeq" = NOISeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot2.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    alpha = 0.50,
                    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                            1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)
                    
venn.plot <- venn.diagram(x=list("DESeq2" = DESeq2.sig.gene, "DESeq" = DESeq.sig.gene, "baySeq" = baySeq.sig.gene, 
                     "NOISeq" = NOISeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot2_without_SAMSeq.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    alpha = 0.50,
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)

## bar chart for number of significant genes

method <- c("EBSeq", "NOISeq", "baySeq", "DESeq", "Voom", "edgeR Exact", "edgeR GLM", "DESeq2", "SAMSeq" )
count <- c(378, 409, 523, 598, 1037, 1235, 1262, 1325, 1775)
names(count) <- method

colors = c("red", "gold", "lightgreen", "violet", "orange", "lightblue", "pink", "cyan4", "coral3") 

par(mfrow=c(1,1), cex=1, mar=c(7, 4, 2, 2)+0.2)
barplot(count, ylim = c(0, 2000), col=colors, ylab = "Number of significant genes", main = "", las = 2) 



                
                    