
aa = read.delim("../../R6186_Karina/result/countStatTable.txt", sep = "\t", header =TRUE) 
#aa = read.delim("../../R7290/nf_results/result/countStatTable.txt", sep = "\t", header =TRUE) 
aa = data.frame(aa)
aa$pct.trimming = aa$Passed.trimming/aa$Total
aa$pct.mapping = aa$Genome.align/aa$Total
aa$pct.miRNAs = aa$Total.reads.in.feature/aa$Total

chiara.samples = aa
philipp.samples = aa

xx = data.frame(rbind(chiara.samples, philipp.samples))
write.table(xx, file = "../results/R7290_contamination/mapping_ratios_Chiara_Philipp_samples.txt", sep = "\t")


annot = read.delim("../../../../scripts/smallRNA_UMI_src/smallRNA_v54/annotation/ensembl/ensemblCel/coarse_seletion_piRNA.gtf", 
                   header = FALSE)
jj = which(annot$V2=="piRNA")

annot.piRNA = annot[jj, c(1:9)]

write.table(annot.piRNA, file = "../../../../scripts/smallRNA_UMI_src/smallRNA_v54/annotation/ensembl/ensemblCel/cel_piRNA.gtf", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


