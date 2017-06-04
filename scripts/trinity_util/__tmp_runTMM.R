library(edgeR)

rnaseqMatrix = read.table("/mnt/lustre_scratch_2012/Tamer/Hus_plant/output/trimmed_RNA_reads/trinity_out_dir/seqclean/bowtie2_align/eXpress_gene_matrix.not_cross_norm.fpkm.tmp", header=T, row.names=1, com='', check.names=F)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
write.table(exp_study$samples, file="/mnt/lustre_scratch_2012/Tamer/Hus_plant/output/trimmed_RNA_reads/trinity_out_dir/seqclean/bowtie2_align/eXpress_gene_matrix.not_cross_norm.fpkm.tmp.TMM_info.txt", quote=F, sep="\t", row.names=F)
