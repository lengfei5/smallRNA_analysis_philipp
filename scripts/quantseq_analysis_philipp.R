##################################################
## Project: Philipp's small RNA seq data with spike-in 
## Script purpose: control data quality and normalize data using spike-ins
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
library("openxlsx")
require('DESeq2')
RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

### data verision and analysis version
version.Data = 'Quantseq_R8043'
version.analysis = paste0("_", version.Data, "_20190719")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE

# spike.concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100 ## the concentration is amol/mug-of-total-RNA

### Directories to save results
design.file = "../exp_design/NGS_Samples_Philipp_mRNAseq_all.xlsx"
dataDir = "../data/"

resDir = paste0("../results/", version.Data, "/")
tabDir =  paste0(resDir, "tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##########################################
# Import Sample information and table of read counts
# mainly manully 
##########################################
if(file.exists(design.file)){
  design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
  design = data.frame(design, stringsAsFactors = FALSE)
  design = design[which(!is.na(design$Sample.ID)), ]
  
  jj = which(colnames(design) == 'Sample.ID')
  design = design[, c(jj, setdiff(c(1:ncol(design)), jj))]
  
  # select columns to keep
  col2keep = c("Sample.ID", "Genetic.Background", "Stage", "Treatment" )
  design = design[, match(col2keep, colnames(design))]
  colnames(design) = c('SampleID', 'strain', 'stage', 'treatment')
  design$condition = NA
  
  ####
  ## manually prepare the design infos
  ####
  #design$strain[grep('Arabidopsis', design$strain)] = "Ath"
  design$strain[grep('H2O', design$strain)] = "H2O.control"
  
  design$stage[grep("cells", design$stage)] = "2cells"
  design$stage[grep('fold', design$stage)] = "2.3.fold"
  design$stage[grep("-", design$stage)] = "none"
  
  design$condition[which(design$strain=="N2")] = "wt"
  design$condition[which(design$strain=="MLC860")] = "pash1.ts"
  design$condition[which(design$strain=="MLC1795")] = "pash1.ts.mirtron"
  design$condition[which(design$strain=="MLC1726")] = "drosha.pash1.aid.RNAi"
  design$condition[which(design$strain=="MLC1729" & design$treatment == "OP50")] = "drosha.pash1.aid.mirtron"
  design$condition[which(design$strain=="MLC1729" & design$treatment == "pash-1 RNAi")] = "drosha.pash1.aid.pash1.RNAi.mirtron"
  design$condition[which(is.na(design$condition))] = 'none'
  design = design[, -which(colnames(design)=="treatment")]
}

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
Dir_umi = paste0(dataDir, "R8043_htseq_counts_BAMs_umi")
Dir_read = paste0(dataDir, "R8043_htseq_counts_BAMs")

source(RNAfunctions)

aa1 <- list.files(path = Dir_umi, pattern = "*out_umiDedup", full.names = TRUE)
aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

aa2 <- list.files(path = Dir_read, pattern = "*.out", full.names = TRUE)
aa2 = merge.countTables.htseq(aa2)
colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))

## compare read counts vs. umi counts
source(RNAfunctions)

pdfname = paste0(resDir, "readCounts_vs_UMI", version.analysis, ".pdf")
pdf(pdfname, width = 10, height = 6)
compare.readCount.UMI(design, aa, normalized = FALSE)

dev.off()

save(design, aa, file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))

######################################
######################################
## Section: spike-in and piRNA normalization
# optionally double check the data quality with sample comparisons
## save the normalized tables
######################################
######################################
Counts.to.Use = "UMI"
QC.for.cpm = FALSE
EDA.with.normalized.table = FALSE
Compare.UMI.vs.readCounts = FALSE

load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

if(Compare.UMI.vs.readCounts){
  pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
  pdf(pdfname, width = 10, height = 8)
  
  compare.readCount.UMI(design, aa, normalized = TRUE)
  
  dev.off()
}

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, special.column = ".readCount")
}else{
  if(Counts.to.Use == "UMI"){
    all = process.countTable(all=aa, design = design, special.column = "UMI")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene

##########################################
# quality control  
##########################################
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))
  
  index.qc = c(1, 4, 3)
  
  source(RNA_QCfunctions)
  
  pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, "_", Counts.to.Use, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, index.qc])
  dev.off()
  
}

##########################################
# calculate scaling factor and normalization
##########################################
if(EDA.with.normalized.table){
  require(DESeq2)
  samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))
  
  raw = ceiling(as.matrix(all[, (samples.sels+1)]))
  raw[which(is.na(raw))] = 0
  rownames(raw) = all$gene
  
  #source(RNAfunctions)
  dds <- DESeqDataSetFromMatrix(raw, 
                                DataFrame(design[samples.sels, ]), 
                                design = ~ condition + stage)
  lowlyExpressed.readCount.threshold = 20
  dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
  dds <- estimateSizeFactors(dds)
  
  fpm = fpm(dds, robust = TRUE)
  
  if(Save.Tables){
    xx = data.frame(fpm, stringsAsFactors = FALSE)
    write.csv(xx, file = paste0(tabDir, "Table_normalized_for_", Counts.to.Use,  version.analysis, ".csv"), 
              row.names = TRUE)
  }
  
}

########################################################
########################################################
# Section : specify pairwise comparisons for different time points and conditions
# 
########################################################
########################################################
kk = which(design$condition == 'none' & design$stage != "none")
design$condition[kk] = design$strain[kk] 

lowlyExpressed.readCount.threshold = 10
require(DESeq2)
source(RNA_QCfunctions)
#index.qc = c(1, 4)

##########################################
# try to merge all analysis for Q1 and Q2
# bacailly all comparisons were done for the same time points
##########################################
compares = list(list("L1", c("MLC1384", "MT17810")),
                list("2.3.fold", c("pash1.ts.mirtron", "drosha.pash1.aid.pash1.RNAi.mirtron", "drosha.pash1.aid.mirtron")),
                list("L1", c("pash1.ts.mirtron", "drosha.pash1.aid.pash1.RNAi.mirtron", "drosha.pash1.aid.mirtron")),
                list("Gastrulation", c("pash1.ts", "pash1.ts.mirtron")), 
                list("Gastrulation", c("drosha.pash1.aid.RNAi", "drosha.pash1.aid.pash1.RNAi.mirtron"))
                )

for(n in 1:length(compares)){
  
  n = 5
  
  stage.sel = compares[[n]][[1]]
  cond.sel = compares[[n]][[2]]
  compName = paste0(c(stage.sel, cond.sel, "N2"), collapse = "_")
  outDir = paste0(resDir, compName)
  
  cat("time to compare -- ", stage.sel, "\n")
  cat("conditions to compare -- ", cond.sel, "\n")
  cat("output Directory -- ", outDir, "\n")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  # select samples for camparisons
  samples.sels = which(design$condition == "wt")
  for(m in 1:length(cond.sel)) samples.sels = c(samples.sels, which(design$condition == cond.sel[m]))
  samples.sels = intersect(which(design$stage == stage.sel), samples.sels)
  
  # check QC
  pdfname = paste0(outDir, "/Data_QC", version.analysis, "_", Counts.to.Use, "_", compName, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, c(1, 4)])
  # start DE analysis
  dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ condition)
  dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
  dds <- estimateSizeFactors(dds)
  
  cpm = fpm(dds, robust = TRUE)
  colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")
  
  dds = estimateDispersions(dds, fitType = "local")
  plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
  abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)
  
  dds = nbinomWaldTest(dds, betaPrior = TRUE)
  resultsNames(dds)
  
  res = cpm
  index.upper = c()
  index.lower = c()
  padj.cutoff = 0.1
  
  for(ii in 1:length(cond.sel)){
    res.ii = results(dds, contrast=c("condition",cond.sel[ii], "wt"))
    summary(res.ii)
    plotMA(res.ii, ylim = c(-2, 2), main = paste0("MA plot -- ", cond.sel[ii], " vs wt"))
    
    ## special cases:  
    # list("Gastrulation", c("pash1.ts", "pash1.ts.mirtron")), 
    # list("Gastrulation", c("drosha.pash1.aid.RNAi", "drosha.pash1.aid.pash1.RNAi.mirtron"))
    if(stage.sel == "Gastrulation" & (cond.sel[ii] == "pash1.ts.mirtron" | cond.sel[ii] == "drosha.pash1.aid.pash1.RNAi.mirtron")){
      upper.ii = which(res.ii$padj > padj.cutoff)
      lower.ii = upper.ii
    }else{
      upper.ii = which(res.ii$log2FoldChange > 0 &  res.ii$padj < padj.cutoff)
      lower.ii = which(res.ii$log2FoldChange < 0 &  res.ii$padj < padj.cutoff)
    }
    
    if(ii == 1) {
      index.upper = c(index.upper, upper.ii)
      index.lower = c(index.lower, lower.ii)
    }else{
      index.upper = intersect(index.upper, upper.ii)
      index.lower = intersect(index.lower, lower.ii)
    }
    
    colnames(res.ii) = paste0(colnames(res.ii), paste0("_", cond.sel[ii], ".vs.N2"))
    res =  data.frame(res, res.ii[, c(1, 2, 5,6)])
  }
  
  dev.off()
  
  write.csv(res, 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  write.csv(res[index.upper, ], 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_Upregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  write.csv(res[index.lower, ], 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_Lowregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  head(res[index.upper, grep("log2Fold|padj", colnames(res))])
  head(res[index.lower, grep("log2Fold|padj", colnames(res))])
  
}

##########################################
# Q3: Compare two time points, gastrulation vs. 2cell stage, separately for
# N2, two mutants and two rescues 
##########################################
stage.sel = c('2cells', 'Gastrulation')
cond.sel = c('wt','pash1.ts', 'drosha.pash1.aid.RNAi', "pash1.ts.mirtron", "drosha.pash1.aid.pash1.RNAi.mirtron")
compName = paste0(c(stage.sel), collapse = "_vs_")
outDir = paste0(resDir, compName)

cat("time to compare -- ", stage.sel, "\n")
cat("conditions to compare -- ", cond.sel, "\n")
cat("output Directory -- ", outDir, "\n")
if(!dir.exists(outDir)) dir.create(outDir)

# select samples for camparisons
samples.sels = c()
for(cc in cond.sel) {
  for(tt in stage.sel)
  samples.sels = c(samples.sels, which(design$condition == cc & design$stage == tt))
}
samples.sels = unique(samples.sels)

# check QC
pdfname = paste0(outDir, "/Data_QC", version.analysis, "_", Counts.to.Use, "_", compName, ".pdf")
pdf(pdfname, width = 12, height = 10)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, c(1, 4, 3)])

#dev.off()
# start DE analysis
dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ condition)
dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
dds <- estimateSizeFactors(dds)

cpm = fpm(dds, robust = TRUE)
colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")

dds = estimateDispersions(dds, fitType = "local")
plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res = cpm
index.upper = c()
index.lower = c()
padj.cutoff = 0.1

for(ii in 1:length(cond.sel)){
  res.ii = results(dds, contrast=c("condition",cond.sel[ii], "wt"))
  summary(res.ii)
  plotMA(res.ii, ylim = c(-2, 2), main = paste0("MA plot -- ", cond.sel[ii], " vs wt"))
  
  ## special cases:  
  # list("Gastrulation", c("pash1.ts", "pash1.ts.mirtron")), 
  # list("Gastrulation", c("drosha.pash1.aid.RNAi", "drosha.pash1.aid.pash1.RNAi.mirtron"))
  if(stage.sel == "Gastrulation" & (cond.sel[ii] == "pash1.ts.mirtron" | cond.sel[ii] == "drosha.pash1.aid.pash1.RNAi.mirtron")){
    upper.ii = which(res.ii$padj > padj.cutoff)
    lower.ii = upper.ii
  }else{
    upper.ii = which(res.ii$log2FoldChange > 0 &  res.ii$padj < padj.cutoff)
    lower.ii = which(res.ii$log2FoldChange < 0 &  res.ii$padj < padj.cutoff)
  }
  
  if(ii == 1) {
    index.upper = c(index.upper, upper.ii)
    index.lower = c(index.lower, lower.ii)
  }else{
    index.upper = intersect(index.upper, upper.ii)
    index.lower = intersect(index.lower, lower.ii)
  }
  
  colnames(res.ii) = paste0(colnames(res.ii), paste0("_", cond.sel[ii], ".vs.N2"))
  res =  data.frame(res, res.ii[, c(1, 2, 5,6)])
}

dev.off()

write.csv(res, 
          file = paste0(outDir, "/NormalizedTable_DEanalysis_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
          row.names = TRUE)

write.csv(res[index.upper, ], 
          file = paste0(outDir, "/NormalizedTable_DEanalysis_Upregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
          row.names = TRUE)
write.csv(res[index.lower, ], 
          file = paste0(outDir, "/NormalizedTable_DEanalysis_Lowregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
          row.names = TRUE)

head(res[index.upper, grep("log2Fold|padj", colnames(res))])
head(res[index.lower, grep("log2Fold|padj", colnames(res))])











