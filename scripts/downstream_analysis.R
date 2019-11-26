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

sels = which(xx$padj_mutant.vs.wt<0.05)
sels1 =  which(xx$padj_mutant.vs.wt<0.05 & xx$log2FoldChange_mutant.vs.wt>0)
sels2 = which(xx$padj_mutant.vs.wt<0.05 & xx$log2FoldChange_mutant.vs.wt<0)

ratios = c((table(xx$chr)/nrow(xx))[7],
           (table(xx$chr[sels])/length(sels))[7],
           (table(xx$chr[sels1])/length(sels1))[7],
           (table(xx$chr[sels2])/length(sels2))[7])
m = length(which(xx$chr=='X'))
n = length(which(xx$chr != "X"))
pvs = c(phyper(q=(length(which(xx$chr[sels] == 'X'))-1), m, n, length(sels), lower.tail = FALSE, log.p = FALSE), 
        phyper(q=(length(which(xx$chr[sels1] == 'X'))-1), m, n, length(sels1), lower.tail = FALSE, log.p = FALSE), 
        phyper(q=(length(which(xx$chr[sels2] == 'X'))-1), m, n, length(sels2), lower.tail = TRUE, log.p = FALSE))
barplot(ratios, names.arg = c('bg', 
                              paste0('dereg, pval = ', signif(pvs[1], d=2)),
                              paste0('upreg, pval = ', signif(pvs[2], d=2)),
                              paste0('downreg, pval = ', signif(pvs[3], d=2))))

##########################################
# PCA plot for samples of interest
##########################################
Counts.to.Use = "UMI"

load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

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
# estimate sample timing 
##########################################



