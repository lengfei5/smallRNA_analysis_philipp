##################################################
## Project: Philipp's small RNA seq data with spike-in 
## Script purpose: control data quality and normalize data using spike-ins
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
library("openxlsx")
require('DESeq2')
miRNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/miRNAseq_functions.R"

### data verision and analysis version
version.Data = 'miRNAs_R7708'
version.analysis = paste0("_", version.Data, "_20190426")

Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE

spike.concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100 ## the concentration is amol/mug-of-total-RNA

### Directories to save results
design.file = "../exp_design/sampleInfo_R7708.xlsx"
dataDir = "../data/"

resDir = "../results/R7708_timeSeries/"
tabDir =  paste0(resDir, "tables/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

##########################################
# Import Sample information and table of read counts
# mainly manully 
##########################################
if(file.exists(design.file)){
  design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
  design = data.frame(design, stringsAsFactors = FALSE)
  jj = which(colnames(design) == 'Sample.ID')
  design = design[, c(jj, setdiff(c(1:ncol(design)), jj))]
  design = design[, c(1:4, 8)]
  colnames(design) = c('SampleID', 'strain', 'stage', 'treatment', 'adaptor')
  design$strain[grep('Arab', design$strain)] = "Ath"
  design$treatment[grep("-", design$treatment)] = 'notreat'
  design$adaptor[grep('Stand', design$adaptor)] = "standard"
  
  #design = design[, -2]
  
}else{
  design = data.frame(
    c(seq(81612, 81615), seq(82086, 82088)), 
    c(rep("cel", 5), "arab", "h2o"),
    c(rep("2cells", 2), rep("L1", 2), "2cells", "None", "None")
  )
  
  colnames(design)[c(1:3)] = c("SampleID", "organism", "stage")
  #design$genotype[grep("WT", design$genotype)] = "WT"
  #design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole_animal_no_promoter"
}

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
aa1 = read.delim(paste0(dataDir, "R7708_result_srbc/countTable.txt"), sep = "\t", header = TRUE)
aa2 = read.delim(paste0(dataDir, "R7708_result_old/countTable.txt"), sep = "\t", header = TRUE)
aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Name", all = TRUE), list(aa1, aa2))

# spike-ins 
spikes1 = read.delim(paste0(dataDir, "R7708_result_srbc/spikeInTable.txt"), sep="\t", header = TRUE, row.names = 1)
spikes2 = read.delim(paste0(dataDir, "R7708_result_old/spikeInTable.txt"), sep="\t", header = TRUE, row.names = 1)
spikes = data.frame(spikes1, spikes2, stringsAsFactors = FALSE)

# stat 
# ignore the stat now, because the nf-pipeline is still yield consistent output
#stat1 = read.delim(paste0(dataDir, "R7708_result_srbc/countStatTable.txt"), sep="\t", header = TRUE)
#stat2 = read.delim(paste0(dataDir, "R7708_result_old/countStatTable.txt"), sep="\t", header = TRUE)
#stat = t(as.matrix(stat))
#spikes = process.countTable(spikes, design)

source(miRNAfunctions)
pdfname = paste0(resDir, "readCounts_vs_UMI", version.analysis, ".pdf")
pdf(pdfname, width = 10, height = 6)
compare.readCounts.umiFr.umiNum(design, aa, spikes)
dev.off()

# start to process table and merge them into one
kk = grep("piRNA_", aa$Name)
if(length(kk)>0){
  piRNAs = aa[kk, ]
  aa = aa[-kk,]
}

source(miRNAfunctions)

Counts.to.Use = "readCounts"

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, select.counts = "Total.count")
}else{
  if(Counts.to.Use == "UMIfr"){
    all = process.countTable(all=aa, design = design, select.counts = "Total.UMIfr.count")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

xx = as.matrix(all[, -1])
xx[which(is.na(xx))] = 0
stat.miRNAs = floor(apply(xx, 2, sum))

if(Counts.to.Use == 'readCounts'){
  piRNAs = process.countTable(all= piRNAs, design = design, select.counts = "GM.count")
}else{
  if(Counts.to.Use == "UMIfr"){
    piRNAs = process.countTable(all= piRNAs, design = design, select.counts = "Total.UMIfr.count")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for piRNAs\n")
  }
}

piRNAs = as.matrix(piRNAs[, -1])
piRNAs[which(is.na(piRNAs))] = 0
stat.piRNAs = floor(apply(piRNAs, 2, sum))
#all = rbind(c('total.piRNAs', stat.piRNAs), all)

spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE)
if(Counts.to.Use == "readCounts"){
  spikes = process.countTable(all=spikes, design = design, select.counts = "Total.spikeIn")
}else{
  if(Counts.to.Use == "UMIfr"){
    spikes = process.countTable(all=spikes, design = design, select.counts = "Total.UMI.spikeIn")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for spike-ins \n")
  }
}

total.spikes = floor(apply(as.matrix(spikes[, -1]), 2, sum))

all = rbind(spikes, all);
stats = rbind(stat.miRNAs, stat.piRNAs, total.spikes)

######################################
######################################
## Section: spike-in and piRNA normalization
# optionally double check the data quality with sample comparisons 
## save the normalized tables
######################################
######################################
read.count = all[, -1];
sel.samples.with.spikeIns = c(1:nrow(design))

raw = floor(as.matrix(read.count[, sel.samples.with.spikeIns]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene
#dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ treatment + stage)
index.spikeIn = grep("spikeIn", rownames(raw))[c(1:8)]


## calculate scaling factor using spike-ins
source("miRNAseq_functions.R")

pdfname = paste0(resDir, "/Spike_in_signals_normalized_DESeq", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 10)
par(mfrow=c(2,2))

res.spike.in = calculate.scaling.factors.using.spikeIns(raw, concentrations = concentrations, index.spikeIn = index.spikeIn, read.threshold = 5)

dev.off()

cpm = res.spike.in$cpm;
res = res.spike.in$normalization.spikeIn
colnames(cpm) = paste0(colnames(cpm), ".cpm")
colnames(res) = paste0(colnames(res), ".amol.per.mugRNA.normBySpikeIns")

ss = apply(raw, 2, sum)
plot(raw[,1]/ss[1]*10^6, cpm[,1], log='xy');abline(0, 1, lwd=2.0, col='red')

#plot(raw[,1]/ss[1]*10^6/norms[1], res[,1], log='xy');abline(0, 1, lwd=2.0, col='red')

### normalization counts using piRNA data
#piRNAs = stat[which(stat$type=="piRNA"), c(-1)]
piRNAs = stat
sizefactors = as.numeric(piRNAs)
cpm.piRNA = raw
for(n in 1:ncol(cpm.piRNA))
{
  cpm.piRNA[,n] = raw[,n]/sizefactors[n]*10^6
}

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

if(Save.Tables){
  colnames(cpm.piRNA) = paste0(colnames(cpm.piRNA), "normBy.piRNA")
  res = data.frame(raw, cpm, res, cpm.piRNA)
  write.csv(res, file = paste0(tabDir, "table_rawCounts_cpm_spikeIn_piRNAs_normalized", version.analysis, ".csv"), 
            row.names = TRUE)
}


########################################################
########################################################
# Section :  check quality
#  this is optional part
########################################################
########################################################
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


######################################
######################################
## Section: check the distribution of normalized data
######################################
######################################
pdfname = paste0(resDir, "Spike_in_signals_normalized_DESeq", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 10)

par(mfrow=c(2,3))
for(n in 1:ncol(cpm))
{
  plot((cpm[kk, n]+10^-6),  concentrations, log='xy', 
       xlab="cpm", ylab="zmol", main=colnames(cpm)[n], cex=1.6);
  abline(h=10^-6, lwd=2.0, col='darkgray')
  points(range((cpm[kk, n]+10^-6)), range((cpm[kk, n]+10^-6))*norms[n], type = 'l', lwd=2.0, col='blue', lty=2)
}

par(mfrow=c(2,3))
for(n in 1:ncol(cpm))
{
  plot((cpm[kk, n]+10^-6),  concentrations, log='', 
       xlab="cpm", ylab="zmol", main=colnames(cpm)[n], cex=1.6);
  abline(h=10^-6, lwd=2.0, col='darkgray')
  points(range((cpm[kk, n]+10^-6)), range((cpm[kk, n]+10^-6))*norms[n], type = 'l', lwd=2.0, col='blue', lty=2)
}

par(mfrow=c(1,1))
xx = log2(res+2^-6)
library(vioplot)
vioplot(xx[,1], xx[,2], xx[,3], xx[,4], xx[,5], xx[,6], names=design$genotype, horizontal = FALSE,
        col="gold")
#title("Violin Plots of Miles Per Gallon")

dev.off()

write.table(res, file = paste0(tabDir, "normalized_signals_using_preliminary_scalling_factors_spike_ins.txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

