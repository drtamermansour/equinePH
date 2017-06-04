library(limma, lib.loc = "~/R/v3.0.1/library")
library(edgeR, lib.loc = "~/R/v3.0.1/library")

files <- list.files(".","RM.PE.quant.counts")
counts <- readDGE(files)  #dim(counts$counts) [1] 43417    24
counts$samples
head(counts$samples)

## low expression threathold
keep<-rowSums(cpm(counts)>0.25) >= 5
counts <- counts[keep, , keep.lib.sizes=FALSE]  #dim(counts$counts) [1] 28900    24  // [1] 28452    22

## normalization
counts <- calcNormFactors(counts)
write.csv(counts$samples, "winter_samples.csv")
write.csv(cpm(counts), "winter_counts.csv")

## grouping factor
counts$samples$group <- factor(c(rep("D", 6), rep("H",6)))
#counts$samples$group <- factor(c(rep("D", 5), rep("H",6)))

#estimate dispersion
y <- estimateCommonDisp(counts)
y <- estimateTagwiseDisp(y)

#run statistical analysis
winterDisVsSummerHealth <- exactTest(y)
write.csv(topTags(winterDisVsSummerHealth, n=100), file="winterDisVsSummerHealth.csv", quote=F)

## prepare for clustering
counts_cpm=cpm(counts)
newCol=substr(colnames(counts),6,nchar(colnames(counts))-9)
colnames(counts_cpm)=newCol
write.table(counts_cpm, file="winter_counts.tab", sep='\t', quote=F, row.names=T)

dir.create("winterDisVsSummerHealth", recursive = TRUE)
Condition <- c(rep("disease", 6), rep("healthy",6))
#Condition <- c(rep("disease", 5), rep("healthy",6))
samples_file <- cbind(Condition , newCol)
write.table(samples_file, file=paste("winterDisVsSummerHealth","samples_file",sep="/"), sep='\t', quote=F, row.names=F, col.names=F)
table.name = paste("eXpress_isoform_matrix.","disease","_vs_","healthy",".edgeR.DE_results",sep="")
write.table(topTags(winterDisVsSummerHealth, n=NULL), file=paste("winterDisVsSummerHealth",table.name,sep="/"), sep='\t', quote=F, row.names=T)

