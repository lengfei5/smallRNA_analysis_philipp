##########################################################################
##########################################################################
# Project:
# Script purpose: some downstream analysis for Philipp's project
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 26 13:00:44 2019
##########################################################################
##########################################################################
version.Data = 'Quantseq_R8043_R8521'
version.analysis = paste0("_", version.Data, "_20190926")

resDir = paste0("../results/", version.Data, "/")
RdataDir = paste0(resDir, "/Rdata/")

RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

##########################################
# Q1: test enrichment of genes in chromosomal X for deregulated genes in mutant.vs.wt at gastrulation stage 
##########################################
stage.sel = 'Gastrulation' # 'Gastrulaiton' or '2cell'
outDir = paste0(resDir, "2cells_vs_Gastrulation_mir35KO/tables_", stage.sel)

res.table = paste0(outDir, "/GeneAll_wt_mutant_rescue_mir35KO_inGastrulation_UMI_2cells_vs_Gastrulation_Quantseq_R8043_R8521_20190926.csv")
xx = read.csv(res.table, header = TRUE, row.names = 1)
jj = grep('mutant.vs.wt', colnames(xx))
xx = xx[,jj]
plot(xx[,1], -log10(xx[, 2]))

## import gene annotation to extract chromosome information
load(file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
mm = match(rownames(xx), annot$Gene.name)
xx = data.frame(xx, chr=annot$Chromosome.scaffold.name[mm], stringsAsFactors = FALSE)
xx = xx[which(xx$chr != 'MtDNA'), ]

kk = which(xx$chr == 'X')
plot(xx[kk,1], -log10(xx[kk, 2]), cex =0.7, xlab = 'log2FC', ylab = '-log10(p-adj)', main = 'genes in chrX')

sels = which(xx$padj_mutant.vs.wt<0.05)
sels1 =  which(xx$padj_mutant.vs.wt<0.05 & xx$log2FoldChange_mutant.vs.wt>0)
sels2 = which(xx$padj_mutant.vs.wt<0.05 & xx$log2FoldChange_mutant.vs.wt<0)

rrs = xx[sels1, ]
boxplot(rrs$log2FoldChange_mutant.vs.wt ~ rrs$chr); abline(h=1, col = 'red')


ratios = c()
pvs = c()
chrs.all = c('I', 'II', 'III', 'IV', 'V', 'X')
counts = table(xx$chr)
for(chr in chrs.all)
{
  ii = which(names(counts) == chr)
  ratios = cbind(ratios, c((table(xx$chr)/nrow(xx))[ii],
             (table(xx$chr[sels])/length(sels))[ii],
             (table(xx$chr[sels1])/length(sels1))[ii],
             (table(xx$chr[sels2])/length(sels2))[ii]))
  
  m = length(which(xx$chr==chr))
  n = length(which(xx$chr != chr))
  pvs = cbind(pvs, c(phyper(q=(length(which(xx$chr[sels] == chr))-1), m, n, length(sels), lower.tail = FALSE, log.p = FALSE), 
          phyper(q=(length(which(xx$chr[sels1] == chr))-1), m, n, length(sels1), lower.tail = FALSE, log.p = FALSE), 
          phyper(q=(length(which(xx$chr[sels2] == chr))-1), m, n, length(sels2), lower.tail = TRUE, log.p = FALSE)))
}
colnames(ratios) = chrs.all
rownames(ratios) = c('bg', 'dereg', 'upreg', 'downreg')
colnames(pvs) = chrs.all
rownames(pvs) = rownames(ratios)[2:4]
barplot(ratios, beside = TRUE, col = c(1:4), legend.text = rownames(ratios), args.legend =list(bty='n') )

write.csv(pvs, paste0(resDir, "downsteam_analysis/enrichment_chrXgenes_pvalues.csv"), col.names = TRUE, row.names = TRUE)

##########################################
# PCA plot for samples of interest
##########################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

Counts.to.Use = "UMI"

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
}else{
  if(Counts.to.Use == "UMI"){
    all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene

require(DESeq2)
require(ggplot2)
source(RNA_QCfunctions)

library("openxlsx")
samples = read.xlsx(paste0(resDir, "downsteam_analysis/PC_samples.xlsx"), sheet = 1, colNames = TRUE)
mm = match(samples$Sample.ID, design$SampleID)
xx = raw[, mm]
design = design[mm, ]
design$Name = samples$Name

xx = as.matrix(xx)
xx[which(is.na(xx))] = 0
dds <- DESeqDataSetFromMatrix(xx, DataFrame(design), design = ~ stage + condition)
lowlyExpressed.readCount.threshold = 10
dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('stage', 'condition') , returnData = TRUE))
pca2save$name = design$Name

ggplot(data=pca2save, aes(PC1, PC2, color=condition, shape=stage, label = name)) + geom_point(size=3) + 
geom_text(hjust = 0.4, nudge_y = 0.5, size=2.5)
#plot(ggp);


##########################################
# estimate sample timing for quant-seq samples
##########################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

Counts.to.Use = "readCounts"

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
}else{
  if(Counts.to.Use == "UMI"){
    all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene


## import the table downloaded from GEO
ss = apply(raw, 2, sum)
test = raw
for(n in 1:ncol(test)) test[,n]  = test[,n]/ss[n]*10^6

test = log2(test + 2^-6)

estimation = rep(NA, ncol(test))
source("/Volumes/groups/cochella/jiwang/Projects/Aleks/scRNAseq_MS_lineage/scripts_analysis/timingEst_functions.R")


if(fastEstimate){
  
  timerGenes.pval = 1; 
  lineageCorrs = NA; 
  loess.span = 0.3;
  lowFilter.threshold.target = 1;
  use = 'lowFilter.target'
  PLOT.test = TRUE
  
  dataDir.Hashimsholy = '/Volumes/groups/cochella/jiwang/Projects/Aleks/scRNAseq_MS_lineage_4save/scRNAseq_MS_lineage/data/Hashimsholy_et_al'
  load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints_tempCorrelation.Rdata"))
  
  if(!is.na(lineageCorrs)){
    kk = which(tcors$AB> lineageCorrs & tcors$MS > lineageCorrs & tcors$E> lineageCorrs & tcors$C>lineageCorrs & timers$pval.box<timerGenes.pval)
    timers = timers[kk, ]
  }
  
  cat('nb of timer genes after filtering : ', length(which(timers$pval.box<timerGenes.pval)), "\n")
  
  saveName = paste0(resDir, "downsteam_analysis/TimingEst_Quantseq_timerGenes_pv.", timerGenes.pval, '_lineageCorr.', lineageCorrs)
  pdfname = paste0(saveName, ".pdf") 
  pdf(pdfname, width = 10, height = 8)
  
  for(kk in c(1:ncol(test))){
    #kk = 10
    estimation[kk] = fast.estimate.timing.with.timer.genes(vec = test[,kk], timers = timers, timepoints = timepoints,
                                                             PLOT.test = PLOT.test,
                                                             timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                                                             lowFilter.threshold.target = lowFilter.threshold.target)
    if(design[kk, 4] == 'wt')   cat(kk, "--", estimation[kk], "--",  design[kk, 3], "-", design[kk, 4], "\n")
  }
  
  dev.off()
  
  design$timingEst = estimation
  write.csv(design, file = paste0(saveName, '.csv'), col.names = TRUE, row.names = FALSE)
  
}




##########################################
# check the introns and exons 
##########################################
##########################################
# estimate sample timing 
##########################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

xx = design[which(design$stage == '2cells'|design$stage == 'Gastrulation'), ]
kk = grep('wt|RNAi|mir35.ko', xx$condition)
xx = xx[kk, ]


