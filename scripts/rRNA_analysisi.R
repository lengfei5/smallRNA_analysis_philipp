##########################################################################
##########################################################################
# Project: Philipp's small RNA 
# Script purpose: compare the rRNA and mRNAs
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Feb 17 15:03:18 2020
##########################################################################
##########################################################################
library("openxlsx")
require('DESeq2')

RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

### data verision and analysis version
version.Data = 'Quantseq_R8043_R8521_rRNA'
version.analysis = paste0("_", version.Data, "_20200217")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE

# spike.concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100 ## the concentration is amol/mug-of-total-RNA

### Directories to save results
design.file = "../exp_design/rRNA_sampleInfos.xlsx"
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
  
  Extract.Emilio.Paula.data = FALSE
  if(Extract.Emilio.Paula.data){
    ## Emilio and Paula's data
    jj = which(design$SampleID >= 100003 | design$strain == "MLC1384" | design$strain == "MT17810")
    design$condition[jj] = design$treatment[jj] 
  }else{
    ## just select Philipp's data
    kk = which(design$SampleID < 100003 & design$strain != "MLC1384" & design$strain != "MT17810")
    design = design[kk, ]
    
  }
  
  ####
  ## manually prepare the design infos
  ####
  #design$strain[grep('Arabidopsis', design$strain)] = "Ath"
  #design$strain[grep('H2O', design$strain)] = "H2O.control"
  #design = design[-which(design$strain == "Ath" | design$strain == "H2O.control"), ]
  
  design$stage[grep("cells", design$stage)] = "2cells"
  design$stage[grep('fold', design$stage)] = "2.3.fold"
  #design$stage[grep("-", design$stage)] = "none"
  
  design$condition[which(design$strain=="N2")] = "wt"
  design$condition[which(design$strain=="MLC860")] = "pash1.ts"
  design$condition[which(design$strain=="MLC1795")] = "pash1.ts.mirtron"
  design$condition[which(design$strain=="MLC1726")] = "drosha.pash1.aid.RNAi"
  design$condition[which(design$strain=="MLC1729" & design$treatment == "OP50")] = "drosha.pash1.aid.mirtron"
  design$condition[which(design$strain=="MLC1729" & design$treatment == "pash-1 RNAi")] = "drosha.pash1.aid.pash1.RNAi.mirtron"
  
  design$condition[which(design$strain == "MT14533" & design$treatment == "20 degree")] = 'mir35.ko.20degree'
  design$condition[which(design$strain == "MT14533" & design$treatment == "25 degree")] = 'mir35.ko.25degree'
  #design$condition[which(is.na(design$condition))] = 'none'
  
  design = design[, -which(colnames(design)=="treatment")]
  
}

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
Dir_umi = paste0(dataDir, "htseq_counts_BAMs_umi_rRNA")
#Dir_read = paste0(dataDir, "R8043_R8521_htseq_counts_BAMs")

source(RNAfunctions)

aa1 <- list.files(path = Dir_umi, pattern = "*out_umiDedup.txt", full.names = TRUE)

aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

#aa2 <- list.files(path = Dir_read, pattern = "*.out", full.names = TRUE)
#aa2 = merge.countTables.htseq(aa2)
#colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

#aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))
aa = aa1
## compare read counts vs. umi counts
# source(RNAfunctions)
# Compare.UMI.vs.readCounts = FALSE
# if(Compare.UMI.vs.readCounts){
#   pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
#   pdf(pdfname, width = 10, height = 8)
#   
#   compare.readCount.UMI(design, aa, normalized = TRUE)
#   
#   dev.off()
# }

#save(design, aa, file = paste0(RdataDir, 'Design_Raw_readCounts_UMI_All_incl_Paula_Emilio', version.analysis, '.Rdata'))
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
#EDA.with.normalized.table = FALSE
#Add.miRNA.Targets = TRUE

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

# all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene

ss = apply(raw, 1, sum)
raw = raw[ss>0, ]

annot = read.csv(file = "/Volumes/groups/cochella/jiwang//annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)

ss = apply(raw, 2, sum)
for(t in c('protein_coding', 'rRNA'))
{
  sels = annot$Gene.name[which(annot$Gene.type==t)]
  kk = match(as.character(sels), rownames(raw))
  kk = unique(kk[which(!is.na(kk))])
  
  xx = raw[kk, ]
  ss = rbind(ss, apply(xx, 2, sum))
  #ss = rbind(ss, apply(xx, 2, sum)/ss[1, ])
  #rr = apply(xx, 2, sum)/apply(raw, 2, sum)
}

rownames(ss) = c('total.reads', 'protein.coding.reads', 'rRNA.reads')

#ss = rbind(ss, ss[,])
write.csv(ss, file = paste0(tabDir, 'stat_mappedReads_proteinCoding_rRNA.csv'), col.names = TRUE, row.names = TRUE)


