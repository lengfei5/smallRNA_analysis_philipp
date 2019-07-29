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

pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
pdf(pdfname, width = 10, height = 6)
compare.readCount.UMI.normalized(design, aa, normalized = TRUE)

dev.off()

save(design, aa, file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))

######################################
######################################
## Section: spike-in and piRNA normalization
# optionally double check the data quality with sample comparisons
## save the normalized tables
######################################
######################################
Counts.to.Use = "readCounts"
QC.for.cpm = TRUE


load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)

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
# xx = as.matrix(all[, -1])
# xx[which(is.na(xx))] = 0
samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))

raw = ceiling(as.matrix(all[, (samples.sels+1)]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene
#raw = raw[which()] 
#dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ treatment + stage)

##########################################
# quality control  
##########################################
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  index.qc = c(1, 2, 3)
  
  source(RNA_QCfunctions)
  
  pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, "_", Counts.to.Use, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[samples.sels, index.qc])
  dev.off()
  
}

##########################################
# calculate scaling factor and normalization
##########################################
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


########################################################
########################################################
# Section :  check quality
#  this is optional part
########################################################
########################################################
if(check.quality.by.sample.comparisons){
  pdfname = paste0(resDir, "/Spike_in_vs_piRNAs", version.analysis, ".pdf")
  pdf(pdfname, width = 10, height = 6)
  
  par(mfrow=c(2,2))
  for(n in 1:2)
  {
    plot(res[, c((n*2-1), (2*n))], log='xy', main = 'spikeIns.norm', xlab = colnames(res)[(2*n-1)],
         ylab = colnames(res)[2*n], cex=0.7)
    abline(0, 1, lwd=2.0, col = 'red')
    
    plot(cpm.piRNA[, c((n*2-1), (2*n))], log='xy', main = 'piRNA.norm', xlab = colnames(cpm.piRNA)[(2*n-1)],
         ylab = colnames(cpm.piRNA)[2*n], cex=0.7)
    abline(0, 1, lwd=2.0, col = 'red')
  }
  
  par(mfrow=c(1,1))
  plot(res.spike.in$norms4DESeq2, piRNAs, log='xy', main = "spikeIns norm vs piRNAs norm")
  text(res.spike.in$norms4DESeq2, piRNAs, labels = colnames(raw), offset = 0., pos = 1)
  #plot(raw[,1]/prod(ss)^(1/length(ss))*10^6/si, res[,1], log='xy');abline(0, 1, lwd=2.0, col='red')
  
  matplot(concentrations, cpm.piRNA[index.spikeIn, ], type = 'b', log='xy', col = c(1:ncol(cpm.piRNA)),
          xlab = "spike concentration (amol/mug total RNA)", ylab='normalized by piRNAs', lwd=1.5, pch = 16)
  legend("topleft", inset=0.01, legend=colnames(cpm.piRNA), col=c(1:ncol(cpm.piRNA)),
         pch=16,  bg= ("white"), horiz=F)
  
  rr = median(c(cpm.piRNA[index.spikeIn, 3]/cpm.piRNA[index.spikeIn, 4], 
                cpm.piRNA[index.spikeIn, 3]/cpm.piRNA[index.spikeIn, 1],
                cpm.piRNA[index.spikeIn, 3]/cpm.piRNA[index.spikeIn, 3]))
  
  res[,3] = res[,3]*rr
  
  par(mfrow=c(2,2))
  for(n in 1:2)
  {
    plot(res[, c((n*2-1), (2*n))], log='xy', main = 'spikeIns.norm.new', xlab = colnames(res)[(2*n-1)],
         ylab = colnames(res)[2*n], cex=0.7)
    abline(0, 1, lwd=2.0, col = 'red')
    
    plot(cpm.piRNA[, c((n*2-1), (2*n))], log='xy', main = 'piRNA.norm', xlab = colnames(cpm.piRNA)[(2*n-1)],
         ylab = colnames(cpm.piRNA)[2*n], cex=0.7)
    abline(0, 1, lwd=2.0, col = 'red')
  }
  
  dev.off()
} 


read.count = all[, -1];

kk = c(1:nrow(design))

design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
raw = as.matrix(read.count[,kk])
raw[which(is.na(raw))] = 0
### start enrichment analysis 
raw = floor(raw)
rownames(raw) = all$gene

index.qc = c(1, 3)

#norms = as.numeric(res.spike.in$norms4DESeq2)
norms = as.numeric(piRNAs)/median(as.numeric(piRNAs))

source("RNAseq_Quality_Controls.R")
pdfname = paste0(resDir, "Data_qulity_assessment", version.analysis, ".pdf")

pdf(pdfname, width = 12, height = 10)
Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix[, index.qc])
dev.off()


