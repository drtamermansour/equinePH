library(limma, lib.loc = "~/R/v3.0.1/library")
library(edgeR, lib.loc = "~/R/v3.0.1/library")
#library(Biobase)
#source("/mnt/ls15/scratch/users/mansourt/Tamer/MFM/scripts/trinity_util/R/rnaseq_plot_funcs.R")

files <- list.files(".","counts")
counts <- readDGE(files)  #dim(counts$counts) [1] 43417    24
counts$samples
head(counts$samples)

## low expression threathold
keep<-rowSums(cpm(counts)>0.25) >= 12
counts <- counts[keep, , keep.lib.sizes=FALSE]  #dim(counts$counts) [1] 28900    24  // [1] 28452    22

## normalization
counts <- calcNormFactors(counts)
write.csv(counts$samples, "samples.csv")
write.csv(cpm(counts), "counts.csv")

#design matrix Diseased State (D=diseased, H=healthy) Season (W=winter, S=summer)
#Condition <- c(rep("D", 12), rep("H",12))
#Subject <- c("1","1","2","2","3","3","4","4","5","5","6","6","1","1","2","2","3","3","4","4","5","5","6","6")
#Season <- c("W","S","W","S","W","S","W","S","W","S","W","S","W","S","W","S","W","S","S","W","S","W","S","W")
Condition <- c(rep("D", 10), rep("H",12))
Subject <- c("1","1","2","2","3","3","4","4","5","5","1","1","2","2","3","3","4","4","5","5","5","5")
Season <- c("W","S","W","S","W","S","W","S","W","S","W","S","W","S","W","S","S","W","S","W","S","W")

Targets <- cbind(Subject , Condition , Season)
Targets<-as.data.frame(Targets)

#create factors
Condition <- factor(Targets$Condition, levels=c("H","D"))
Season <- factor(Targets$Season, levels=c("W","S"))

design <- model.matrix(~Condition + Condition:Subject + Condition:Season)
colnames(design)

#estimate dispersion
y <- estimateGLMCommonDisp(counts,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

#fit to a linear model
fit <- glmFit(y, design)

#genes that respond differently in disease state (comparison1)
# represents the log-fold change in diseased winter of patient 1 over healthy winter of patient 7 (useless)
#disVsHealthy <- glmLRT(fit, coef="ConditionD")

#gene changes in disease state in summer vs winter (comparison 2)
# represents the log-fold change in summer diseased over winter diseased
# i.e. effect of summer in the diseased patients
effSummerInDis <- glmLRT(fit, coef="ConditionD:SeasonS")
write.csv(topTags(effSummerInDis, n=100), file="effSummerInDis.csv", quote=F)

#genes chanegs in healthy state, summer vs winter (comparison 3)
# represents the log-fold change in summer healthy over winter healthy
# i.e. effect of summer in the healthy patients
effSummerInHealthy <- glmLRT(fit, coef="ConditionH:SeasonS")
write.csv(topTags(effSummerInHealthy, n=100), file="effSummerInHealthy.csv", quote=F)

#genes that are affected by season depending on disease state (?) (comparison 4)
#lrt <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1))

#genes that respond to summer in healthy and disease (comparison5 is just a combined table for 2 & 3)
#InDisAndHealth_SummerVsWinter <- glmLRT(fit, coef=13:14)

#comparison 4 reversed (comparison 6)
#tests for differences in the summer/winter log-fold change between diseased and healthy patients.
#This allows you to determine whether the effect of summer on transcription is itself affected by the disease.
# i.e. the effect of disease on summer/winter log-fold change
# the best approximation for the overall effect of the disease
#effDisOnSeasonalChange <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1))
effDisOnSeasonalChange <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,1,-1))
write.csv(topTags(effDisOnSeasonalChange, n=100), file="effDisOnSeasonalChange.csv", quote=F)

#genes respond in a overall affect to summer (comparison 7)
#represents the log-fold change in avearge summer change over winter (in diseased and healthy)
#i.e. overall effect of summer in healthy and diseased patients
#SummerVsWinter <- glmLRT(fit2, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2))


## prepare for clustering
counts_cpm=cpm(counts)
newCol=substr(colnames(counts),6,nchar(colnames(counts))-9)
colnames(counts_cpm)=newCol
write.table(counts_cpm, file="counts.tab", sep='\t', quote=F, row.names=T)

dir.create("effDisOnSeasonalChange", recursive = TRUE)
#Condition <- c(rep("disease", 12), rep("healthy",12))
Condition <- c(rep("disease", 10), rep("healthy",12))
samples_file <- cbind(Condition , newCol)
write.table(samples_file, file=paste("effDisOnSeasonalChange","samples_file",sep="/"), sep='\t', quote=F, row.names=F, col.names=F)
table.name = paste("eXpress_isoform_matrix.","disease","_vs_","healthy",".edgeR.DE_results",sep="")
write.table(topTags(effDisOnSeasonalChange, n=NULL), file=paste("effDisOnSeasonalChange",table.name,sep="/"), sep='\t', quote=F, row.names=T)

dir.create("effSummerInDis", recursive = TRUE)
#Condition <- c(rep(c("nonExc","Exacerb"), 6), rep("nonExc",12))
Condition <- c(rep(c("nonExc","Exacerb"), 5), rep("nonExc",12))
samples_file <- cbind(Condition , newCol)
write.table(samples_file, file=paste("effSummerInDis","samples_file",sep="/"), sep='\t', quote=F, row.names=F, col.names=F)
table.name = paste("eXpress_isoform_matrix.","Exacerb","_vs_","nonExc",".edgeR.DE_results",sep="")
write.table(topTags(effSummerInDis, n=NULL), file=paste("effSummerInDis",table.name,sep="/"), sep='\t', quote=F, row.names=T)

dir.create("effSummerInHealthy", recursive = TRUE)
#Condition <- c(rep(c("nonExc","Exacerb"), 6), rep("nonExc",12))
Condition <- c(rep(c("nonExc","Exacerb"), 5), rep("nonExc",12))
samples_file <- cbind(Condition , newCol)
write.table(samples_file, file=paste("effSummerInHealthy","samples_file",sep="/"), sep='\t', quote=F, row.names=F, col.names=F)
table.name = paste("eXpress_isoform_matrix.","Exacerb","_vs_","nonExc",".edgeR.DE_results",sep="")
write.table(topTags(effSummerInHealthy, n=NULL), file=paste("effSummerInHealthy",table.name,sep="/"), sep='\t', quote=F, row.names=T)

## another heatmap approach
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("Blue", "white", "Red"), bias = 25)(n = 299)

#pdf("heatmap.pdf")
#heatmap.2(cpm(counts),    # data matrix

keep<-rowSums(cpm(counts)>1) >= 6
counts2 <- counts[keep, ]      ## 26440
pdf("heatmap2.pdf")
heatmap.2(cpm(counts2),    # data matrix
#cellnote = mat_data,  # same data set for cell labels
main = "heatmap", # heat map title
#notecex=0.4,
cexRow=0.6,
cexCol=0.8,
scale="none",
#notecol="black",      # change font color of cell labels to black
density.info="none",  # turns off density plot inside color legend
trace="none",         # turns off trace lines inside the heat map
margins =c(12,12),     # widens margins around plot
col=my_palette,       # use on color palette defined earlier
#breaks=col_breaks,    # enable color transition at specified limits
#dendrogram="row",     # only draw a row dendrogram
Colv=TRUE)            # turn off column clustering
dev.off()               # close the PNG device



