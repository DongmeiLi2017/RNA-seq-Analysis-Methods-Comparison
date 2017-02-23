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

bottomly_c <- bottomly_count[, sample(2:11, 3)]   ##9 , 5, 11
bottomly_t <- bottomly_count[, sample(12:22, 3)]   ## 22, 12, 16 
#> dim(bottomly_count)
#[1] 36536    22

bottomly_count2 <- cbind(bottomly_count[, 1], bottomly_c, bottomly_t)
#> dim(bottomly_count2)
#[1] 36536     7


keep <- rowSums(bottomly_count2[, -1])>=10
nkeep <- sum(keep)
#> nkeep
#[1] 10815
bottomly_countnew <- bottomly_count2[keep,]
label <- c(rep(0, 3), rep(1, 3))
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
#[1]  143
edgeR.bottomly <- cbind(bottomly_countnew, edgeR.test$table, edgeR.adjp)
edgeR.bottomly.select <- edgeR.bottomly[which(edgeR.adjp <=0.05), ]

edgeR.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(edgeR.adjp <=0.05)])

## edgeR GLM Test 3.12.0

edgeR.glm <- estimateGLMTrendedDisp(edgeR.dgelist, design)
edgeR.glm <- estimateGLMTagwiseDisp(edgeR.glm, design)
edgeR.glmf <- glmFit(edgeR.glm, design)
edgeR.glmde <- glmLRT(edgeR.glmf)
edgeR.glmpvalues <- edgeR.glmde$table$PValue
edgeR.glmadjp <- p.adjust(edgeR.glmpvalues, "BH")
sum(edgeR.glmadjp <=0.05)
#[1]  179
edgeR.glm.bottomly <- cbind(bottomly_countnew, edgeR.glmde$table, edgeR.glmadjp)
edgeR.glm.bottomly.select <- edgeR.glm.bottomly[which(edgeR.glmadjp <=0.05), ]

edgeR.glm.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(edgeR.glmadjp <=0.05)])

## DESeq  1.22.1

bottomly.norm <- round(edgeR.dgelist$counts, 0)
DESeq.cds <- newCountDataSet(countData = bottomly.norm, conditions = factor(label))
DESeq.disp <- estimateSizeFactors(DESeq.cds)
DESeq.disp <- estimateDispersions(DESeq.disp, sharingMode = "maximum", method = "pooled", fitType = "local")
DESeq.test <- nbinomTest(DESeq.disp, "0", "1")
DESeq.pval <- DESeq.test$pval
DESeq.adjp <- p.adjust(DESeq.pval, "BH")
sum(DESeq.adjp <=0.05)
#[1] 81
DESeq.bottomly <- cbind(bottomly_countnew, DESeq.test, DESeq.adjp)
DESeq.bottomly.select <- DESeq.bottomly[which(DESeq.adjp <=0.05), ]

DESeq.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(DESeq.adjp <= 0.05)])

## DESeq2  1.10.1

group.con <- data.frame(cbind(c(rep("C57BL", 3), rep("DBA", 3)), c(rep("single-read", 6))))
colnames(group.con) <- c("condition", "type")
rownames(group.con) <- c(paste("C57BL", 1:3), paste("DBA", 1:3))
colnames(bottomly.norm) <- rownames(group.con)
DESeq2.dds <- DESeqDataSetFromMatrix(countData = bottomly.norm, colData = group.con, design = ~condition)
DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
DESeq2.pval <- results(DESeq2.test)$pvalue
DESeq2.adjp <- p.adjust(DESeq2.pval, "BH")
sum(DESeq2.adjp <= 0.05, na.rm=TRUE)
#[1]  100
DESeq2.bottomly <- cbind(bottomly_countnew, results(DESeq2.test), DESeq2.adjp)
DESeq2.bottomly.select <- DESeq2.bottomly[which(DESeq2.adjp <= 0.05), ]

DESeq2.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(DESeq2.adjp <= 0.05)])

## baySeq 2.4.1

baySeq.CD <- new("countData", data = bottomly.norm, replicates = label, groups =list(NDE = c(rep(1, 6)), DE = label))
libsizes(baySeq.CD) <-getLibsizes(baySeq.CD, estimationType = "edgeR")
baySeq.CD <- getPriors.NB(baySeq.CD, samplesize = 5000, equalDispersions = TRUE, estimation = "QL", cl = NULL)
baySeq.CD <- getLikelihoods(baySeq.CD, cl = NULL, verbose=FALSE)
baySeq.posteriors.DE <- exp(baySeq.CD@posteriors)[, 2]
baySeq.table <- topCounts(baySeq.CD, group = "DE", FDR = 1)
baySeq.FDR <- baySeq.table$FDR[match(rownames(bottomly.norm), rownames(baySeq.table))]
sum(baySeq.FDR <= 0.05, na.rm=TRUE)
# [1] 55
baySeq.bottomly <- cbind(bottomly_countnew, baySeq.table)
baySeq.bottomly.select <- baySeq.bottomly[which(baySeq.FDR <= 0.05), ]

baySeq.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(baySeq.FDR <= 0.05)])

## EBSeq  1.10.0

sizes <- MedianNorm(bottomly.norm)
EBSeq.out <- EBTest(Data = bottomly.norm, Conditions = factor(label), sizeFactors = sizes, maxround = 5)
EBSeq.result <- GetDEResults(EBSeq.out, FDR = 0.05)
EBSeq.status <- EBSeq.result$Status
sum(EBSeq.status == "DE", na.rm=TRUE)
#[1] 108

EBSeq.bottomly <- cbind(bottomly_countnew, EBSeq.out$PPDE, EBSeq.out$f0, EBSeq.out$f1, EBSeq.status)
EBSeq.bottomly.select <- EBSeq.bottomly[which(EBSeq.status == "DE"), ]

EBSeq.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(EBSeq.status == "DE")])

## SAMSeq   2.0
label2<-c(rep(1, 3), rep(2, 3))
SAMSeq.test <- SAMseq(bottomly.norm, label2, resp.type = "Two class unpaired", geneid = rownames(bottomly.norm), genenames = rownames(bottomly.norm), nperms = 100, nresamp = 20, fdr.output = 0.05)
SAMSeq.result <- rbind(SAMSeq.test$siggenes.table$genes.up, SAMSeq.test$siggenes.table$genes.lo)
SAMSeq.s <- strsplit(SAMSeq.result[, 1], "[^[:digit:]]")
SAMSeq.solution <- as.numeric(unlist(SAMSeq.s))
SAMSeq.solution <- unique(SAMSeq.solution[!is.na(SAMSeq.solution)])
length(SAMSeq.solution)
#[1] 206
SAMSeq.bottomly.select <- bottomly_count2[SAMSeq.solution, ]

SAMSeq.sig.gene <- SAMSeq.solution

## NOISeq    2.14.1
NOISeq.data <- readData(bottomly.norm, factors = group.con)
NOISeq.result <- noiseqbio(NOISeq.data, norm = "tmm", factor = "condition", lc = 1, filter = 0)
NOISeq.DE <- degenes(NOISeq.result, q = 0.95, M = NULL)
NOISeq.s <- strsplit(rownames(NOISeq.DE), "[^[:digit:]]")
NOISeq.solution <- as.numeric(unlist(NOISeq.s))
NOISeq.solution <- unique(NOISeq.solution[!is.na(NOISeq.solution)])
length(NOISeq.solution)
#[1] 120
#NOISeq.bottomly <- cbind(bottomly_countnew, NOISeq.result)
#NOISeq.bottomly.select <- EBSeq.bottomly[which(EBSeq.status == "DE"), ]

NOISeq.bottomly.select <- bottomly_count2[NOISeq.solution, ]

NOISeq.sig.gene <- NOISeq.solution

##voom + limma  3.26.8

Voom.data <- voom(bottomly.norm, design = design)
Voom.limma <- lmFit(Voom.data, design, na.rm=TRUE)
Voom.ebayes<- eBayes(Voom.limma)
Voom.pval<-Voom.ebayes$p.value[, 2]
Voom.adjp<-p.adjust(Voom.pval, method = "BH")
sum(Voom.adjp <= 0.05, na.rm=TRUE)
#[1] 43

Voom.bottomly <- cbind(bottomly_countnew, Voom.ebayes, Voom.adjp)
Voom.bottomly.select <- Voom.bottomly[which(Voom.adjp <= 0.05), ]

Voom.sig.gene <- as.numeric(rownames(bottomly_countnew)[which(Voom.adjp <= 0.05)])

write.csv(edgeR.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_edgeR_select_n_6.csv")
write.csv(edgeR.glm.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_edgeR_glm_select_n_6.csv")
write.csv(DESeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_DESeq_select_n_6.csv")
write.csv(DESeq2.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_DESeq2_select_n_6.csv")
write.csv(baySeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_baySeq_select_n_6.csv")
write.csv(EBSeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_EBSeq_select_n_6.csv")
write.csv(SAMSeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_SAMSeq_select_n_6.csv")
write.csv(NOISeq.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_NOISeq_select_n_6.csv")
write.csv(Voom.bottomly.select, file = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_Voom_select_n_6.csv")


select.common <- calculate.overlap(x=list("edgeR_exact" = edgeR.sig.gene, "edgeR_glm" = edgeR.glm.sig.gene))

$a1
  [1]    96   261   506   600   889   935   957  1004  1350  1623  1676  2082  2278  2645  2836  2929  3162
 [18]  3380  3484  3488  3990  4125  4277  4337  4341  4424  4508  4890  5547  5560  5577  5749  6279  6391
 [35]  6483  6669  6703  6816  6845  7077  7270  7294  7573  7582  7936  8031  8447  8678  8828  8970  9029
 [52]  9336  9433  9461  9535  9669  9909 10048 10082 10200 10224 10354 10406 10571 10593 10860 11157 11221
 [69] 11368 11573 11601 11638 11741 11885 12096 12202 12460 12466 12481 12608 13024 13347 13637 13770 13794
 [86] 13809 14215 14342 14564 14647 15017 15036 15213 15453 15517 15612 15666 15740 15977 16106 16108 16165
[103] 16227 16232 16247 16726 17131 17649 18111 18402 18691 18720 20235 20417 20597 20754 20942 20952 21105
[120] 21641 21709 21744 21792 21890 21943 22112 22119 22285 22558 22678 22752 22769 22916 22971 23511 25615
[137] 25616 25677 25740 25853 26154 26431 26499

$a2
  [1]    13    96   261   506   600   889   935   957  1004  1288  1350  1623  1676  1683  2082  2278  2645
 [18]  2836  2929  3162  3380  3484  3488  3990  4125  4277  4337  4341  4424  4508  4546  4841  4890  5547
 [35]  5560  5577  5581  5749  6126  6279  6391  6483  6669  6703  6816  6845  7077  7270  7294  7573  7582
 [52]  7599  7793  7936  8031  8372  8447  8449  8595  8678  8814  8828  8970  8972  9029  9336  9433  9461
 [69]  9535  9669  9770  9909 10047 10048 10082 10200 10224 10354 10406 10571 10593 10860 11157 11221 11368
 [86] 11573 11601 11638 11741 11885 12096 12187 12202 12212 12403 12460 12466 12481 12608 12676 12919 13024
[103] 13347 13637 13658 13770 13794 13809 14215 14342 14386 14564 14647 14858 14867 15017 15036 15213 15453
[120] 15517 15612 15666 15740 15892 15977 16106 16108 16165 16227 16232 16247 16589 16726 17131 17649 17719
[137] 18111 18402 18691 18720 20235 20417 20597 20754 20912 20942 20952 21105 21404 21497 21641 21709 21744
[154] 21792 21890 21943 22112 22119 22132 22226 22285 22558 22678 22752 22769 22795 22916 22971 23511 25615
[171] 25616 25677 25740 25853 25995 26154 26358 26431 26499

$a3
  [1]    96   261   506   600   889   935   957  1004  1350  1623  1676  2082  2278  2645  2836  2929  3162
 [18]  3380  3484  3488  3990  4125  4277  4337  4341  4424  4508  4890  5547  5560  5577  5749  6279  6391
 [35]  6483  6669  6703  6816  6845  7077  7270  7294  7573  7582  7936  8031  8447  8678  8828  8970  9029
 [52]  9336  9433  9461  9535  9669  9909 10048 10082 10200 10224 10354 10406 10571 10593 10860 11157 11221
 [69] 11368 11573 11601 11638 11741 11885 12096 12202 12460 12466 12481 12608 13024 13347 13637 13770 13794
 [86] 13809 14215 14342 14564 14647 15017 15036 15213 15453 15517 15612 15666 15740 15977 16106 16108 16165
[103] 16227 16232 16247 16726 17131 17649 18111 18402 18691 18720 20235 20417 20597 20754 20942 20952 21105
[120] 21641 21709 21744 21792 21890 21943 22112 22119 22285 22558 22678 22752 22769 22916 22971 23511 25615
[137] 25616 25677 25740 25853 26154 26431 26499

                                   
                                   
select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "edgeR_glm" = which(edgeR.glmadjp <=0.05), "DESeq" = which(DESeq.adjp <= 0.05),
                                   "DESeq2" = which(DESeq2.adjp <= 0.05), "baySeq" = which(baySeq.FDR <= 0.05), "EBSeq"= which(EBSeq.status == "DE"),
                                   "SAMSeq" = SAMSeq.solution, "NOISeq" = NOISeq.solution, "Voom" = which(Voom.adjp <= 0.05)))

select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "edgeR_glm" = which(edgeR.glmadjp <=0.05), "DESeq" = which(DESeq.adjp <= 0.05),
                                   "DESeq2" = which(DESeq2.adjp <= 0.05), "baySeq" = which(baySeq.FDR <= 0.05), "EBSeq"= which(EBSeq.status == "DE"),
                                   "SAMSeq" = SAMSeq.solution,  "Voom" = which(Voom.adjp <= 0.05)))

select.common <- calculate.overlap(x=list("edgeR_exact" = which(edgeR.adjp <=0.05), "Voom" = which(Voom.adjp <=0.05), "DESeq2" = which(DESeq2.adjp <= 0.05),
                                          "EBSeq"= which(EBSeq.status == "DE"), "SAMSeq" = SAMSeq.solution))
                                          
select.common <- calculate.overlap(x=list("edgeR_glm" = which(edgeR.glmadjp <=0.05), "edgeR_glm" = which(edgeR.glmadjp <=0.05)))                                         

## venn diagram plot
venn.plot <- venn.diagram(x=list("edgeR Exact" = edgeR.sig.gene, "Voom" = Voom.sig.gene, "DESeq2" = DESeq2.sig.gene,
                    "EBSeq"= EBSeq.sig.gene, "SAMSeq" = SAMSeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot_n_6.png",
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
                    "EBSeq"= EBSeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot_without_SAMSeq_n_6.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    alpha = 0.50,
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)
                    

venn.plot <- venn.diagram(x=list("DESeq2" = DESeq2.sig.gene, "DESeq" = DESeq.sig.gene, "baySeq" = baySeq.sig.gene, 
                    "SAMSeq" = SAMSeq.sig.gene, "NOISeq" = NOISeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot2_n_6.png",
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
                     "NOISeq" = NOISeq.sig.gene), filename = "D:/Dongmei/Research/NGS/RNASeq_examples/bottomly_plot2_without_SAMSeq_n_6.png",
                    col = "black",
                    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    alpha = 0.50,
                    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                    cat.cex = 0.8,
                    cat.fontface = "bold",
                    margin = 0.05)

## bar chart for number of significant genes

method <- c("EBSeq", "NOISeq", "baySeq", "DESeq", "Voom", "edgeR Exact", "edgeR GLM", "DESeq2", "SAMSeq" )
count2 <- c(108, 120, 55, 81, 43, 143, 179, 100, 206)
names(count2) <- method

colors = c("red", "gold", "lightgreen", "violet", "orange", "lightblue", "pink", "cyan4", "coral3") 

par(mfrow=c(1,1), cex=1, mar=c(7, 4, 2, 2)+0.2)
barplot(count2, ylim = c(0, 400), col=colors, ylab = "Number of significant genes", main = "", las = 2) 



                
                    