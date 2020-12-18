##########################################################################
##########################################################################
# Project: Philipp's small RNA-seq 
# Script purpose: generate bigwig files from aligned bam file of small RNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Dec 18 15:17:23 2020
##########################################################################
##########################################################################
rm(list=ls())

library(GenomicAlignments)
library(rtracklayer)

## data verision and analysis version
version.Data = 'R8024_R7846_R7708_all'
version.analysis = paste0("_", version.Data, "_20201218")
resDir = paste0("../results/", version.Data, "/bigWigs/")
#tabDir =  paste0(resDir, "tables/")
#RdataDir = paste0(resDir, "Rdata/")
OutDir = resDir

if(!dir.exists(OutDir)) dir.create(OutDir)

baseDir = '/Volumes/groups/cochella/jiwang/Projects/Philipp/R8024_R7846_R7708_smRNAseq/'

#InputDir = paste0(baseDir, '/Bams_atac_chip_all')
bamlist = list.files(path = paste0(baseDir, 'nf_processing_srbc_v4_makeBW/aligned_bams'), pattern = "*.bam$", full.names = TRUE)
bamlist = unique(c(bamlist, 
            list.files(path = paste0(baseDir, 'nf_processing_old_v4_makeBW/aligned_bams'), pattern = "*.bam$", full.names = TRUE)))

load('/Volumes/groups/cochella/jiwang/Projects/Philipp/smallRNA_analysis_philipp/results/R8024_R7846_R7708_all/Design_Raw_readCounts_UMIfr_miRNAs_R8024_R7846_R7708_20200805.Rdata')

Manual.changeFileName = FALSE
Normalized.coverage = TRUE
Logtransform.coverage = FALSE
Pairend = FALSE
strand.specific = TRUE

for(n in c(1:nrow(design.matrix)))
{
  # n = 1
  bam = bamlist[grep(design.matrix$SampleID[n], bamlist)]
  bw.name = paste0(design.matrix$strain[n], '_', design.matrix$stage[n], '_', design.matrix$treatment[n], '_', design.matrix$SampleID[n])
  bw.name = gsub(" ", ".", bw.name)
  bw.name = gsub("-", '.', bw.name)
  
  if(Manual.changeFileName){
    bw.name = gsub("140min", '200min', bw.name)
    bw.name = gsub("60min", '90min', bw.name)
    bw.name = gsub("Aba", 'ABa', bw.name)
    bw.name = gsub("Abp", 'ABp', bw.name)
  }
  cat(n, '\n')
  cat("bam file: ", bam, "\nbw name: ", bw.name, "\n")
  
  if(!file.exists(paste0(OutDir, bw.name))){
    if(Pairend){
      ga = readGAlignmentPairs(bam)
      #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
    }else{
      ga = readGAlignments(bam)
    }
    
    if(strand.specific){
      alignGR <- granges(ga)
      #strand(gaa)
      alignPos <- alignGR[strand(alignGR) == "+"] 
      alignNeg <- alignGR[strand(alignGR) == "-"] 
      strand(alignPos)
      strand(alignNeg)
      posCov <- coverage(alignPos)
      negCov <- coverage(alignNeg)
      
      if(Normalized.coverage){
        ss = length(alignGR)
        posCov = posCov/ss*10^6
        negCov = negCov/ss*10^6
      }
      
      export.bw(posCov, con = paste0(OutDir, bw.name, '_posCov.bw'))
      export.bw(negCov, con = paste0(OutDir, bw.name, '_negCov.bw'))
      
    }else{
      if(Normalized.coverage){
        ss = length(ga)/2
        xx = coverage(granges(ga))/ss*10^6
      }else{
        xx = ga
      }
      #if(Logtransform.coverage) xx = log2(xx+2^-6)
      export.bw(xx, con = paste0(OutDir, bw.name))
    }
   
  }
}
