find.mirName = function(x)
{
  test = unlist(strsplit(as.character(x), '-'));
  return(paste0(test[-c(1,length(test))], collapse = '-'))
} 

find.mirName.dm = function(x)
{
  test = unlist(strsplit(as.character(x), '-'));
  return(paste0(test[-c(which(test=="RM"),length(test))], collapse = '-'))
} 

find.expressed.mature.miRNA = function(dds, cpm.threshold=10, sampletoUse="untreated")
{
  ## sampletoUse can be "untreated", "treated", "all", "untreated.or.treated"
  require('DESeq2')
  # dds = dds_all
  #index.test = c(294, 332)
  ## filter lowly expressed miRNAs using cpm=10 and select the mature arms always baesed on the untreated samples
  cpm = fpm(dds, robust = FALSE)
  colnames(cpm) = paste0("diluation_", colnames(cpm), '.cpm')
  ggs = sapply(rownames(cpm), find.mirName.dm)
  cpm.untreated = cpm[, which(dds$treatment=="untreated")]
  
  #cpm.treated = cpm[, which(dds$treatment=="treated")]
  if(length(which(dds$treatment=="untreated"))>1) {mean.untreated = apply(cpm[, which(dds$treatment=="untreated")], 1, mean);
  }else {mean.untreated = cpm[, which(dds$treatment=="untreated")];}
  
  if(length(which(dds$treatment=="treated"))>1) {mean.treated = apply(cpm[, which(dds$treatment=="treated")], 1, mean);
  }else {mean.treated = cpm[, which(dds$treatment=="treated")];}
  
  if(sampletoUse=="untreated") gene.expressed = unique(ggs[which(mean.untreated>cpm.threshold)])
  if(sampletoUse=="treated") gene.expressed = unique(ggs[which(mean.treated>cpm.threshold)])
  if(sampletoUse=="all"){
    mean.all = apply(as.matrix(cbind(mean.treated, mean.untreated)), 1, mean);
    gene.expressed = unique(ggs[which(mean.all>cpm.threshold)])
  } 
  if(sampletoUse=="untreated.or.treated") gene.expressed = unique(ggs[which(mean.treated>cpm.threshold|mean.untreated>cpm.threshold)])
  
  mirna.expressed = !is.na(match(ggs, gene.expressed))
  sels = data.frame(miRNA=rownames(cpm), gene=ggs, expressed=mirna.expressed, stringsAsFactors = FALSE)
  
  ggs.uniq = unique(sels$gene)
  cat("nb of expressed genes --", length(gene.expressed), "-- miRNAs", sum(sels$expressed), "\n")
  
  for(n in 1:length(ggs.uniq))
  {
    #n = 1
    jj = which(sels$gene==ggs.uniq[n])
    if(length(jj)>1){
      if(length(which(dds$treatment=="untreated"))>1){index.max = apply(cpm.untreated[jj,], 2, which.max);
      }else{index.max = which.max(cpm.untreated[jj])}
      nb.first.max = length(which(index.max==1))
      nb.sec.max = length(which(index.max==2))
      #cat(n, ": ", as.character(ggs.uniq[n]), "--",  nb.first.max, "--", nb.sec.max, "\n")
      if(nb.first.max>nb.sec.max){ sels$miRNA[jj[1]] = sels$gene[jj[1]];
      }else{sels$miRNA[jj[2]] = sels$gene[jj[2]];}
    }else{
      sels$miRNA[jj] = sels$gene[jj];
    }
  }
  
  return(list(miRNA=sels$miRNA, expressed=sels$expressed, cpm=cpm))
}

Compare.three.Controls.Ovary.dm = function(read.count, design.matrix)
{
  #read.count=read.count[, kk]; design.matrix = design.matrix[, index.qc]
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
  library("ggplot2")
  #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  # kk = grep('Ab', colnames(bb))
  if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
  }else{cc = design.matrix[, -1]}
  #o1 = order(cc)
  #read.count = read.count[o1,]
  #cc = cc[o1]
  raw = as.matrix(read.count)
  #xx = raw
  dim(raw)
  raw[which(is.na(raw))] = 0
  xx = raw;
  
  countData = ceiling(raw)
  conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
  eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
  #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
  #dds <- dds[ rowSums(counts(dds)) > 10, ]
  dds <- estimateSizeFactors(dds)
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  
  ## filter lowly expressed miRNAs
  sels = find.expressed.mature.miRNA(dds)
  rownames(dds) = sels$miRNA;
  dds = dds[sels$expressed, ];
  dds <- estimateSizeFactors(dds)
  #fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  cpm = fpm(dds, robust=TRUE)
  cpm = log2(cpm+0.25)
  
  pairs(cpm[, grep('_untreated', colnames(cpm))], lower.panel=NULL, upper.panel=panel.fitting)
  pairs(cpm[, grep('_treated', colnames(cpm))], lower.panel=NULL, upper.panel=panel.fitting)
  
}


mixture.gaussian.fitting = function(x, method='Normal', k=2, mu.init=NULL, sigma.init=NULL, lambda.init=c(0.3, 0.7), m.constr=NULL, sigma.constr=NULL)
{
    print("to finish")
}

Filtering.Expressed.Genes = function(countData, conds, Use.only.notreated=FALSE, posterior.cutoff=0.5, fraction.detected.samples=0.5)
{
  require('edgeR')
  library(mixtools)
  #cat(design.matrix)
  y0 = countData;
  group = conds
  #group <- apply(design.matrix[, -1], 1, paste0, collapse = "_")
  if(Use.only.notreated){
    jj = which(design.matrix$condition=="notreated");
    y0 = y0[,jj]; group = group[jj];
  }
  y <- DGEList(counts=y0, group=group)
  #y1 <- calcNormFactors(y, method=c('TMM'))
  #y2 <- calcNormFactors(y, method=c('RLE'))
  #y3 <- calcNormFactors(y, method=c('none'))
  y = calcNormFactors(y, method=c("upperquartile"), p=0.5)
  cpm = cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
  
  expressed = matrix(NA, nrow=nrow(cpm), ncol= ncol(cpm));
  par(mfrow=c(2,2))
  #unique.group= unique(group)
  for(n in 1:ncol(cpm)){
    #cc = group[n]
    cat(group[n], '\n');
    data = as.numeric(unlist(cpm[, n]));
    #jj2 = index.outliers(data);
    jj1 = which(y0[, n]>=1);
    jj2 = setdiff(c(1:nrow(expressed)), jj1);
    #samples = design.matrix[which(group==cc), 1]
    #index = c()
    #for(s in samples){index = c(index, which(colnames(cpm)==s));}
    #index = unique(index);
    #if(length(index)>1){data = as.numeric(unlist(apply(cpm[, index], 1, mean)))
    #}else{data = as.numeric(unlist(cpm[, index]));}
    fit = normalmixEM(data[jj1], lambda = c(0.3, 0.7), mu=c(0, 10), sigma = NULL, mean.constr=NULL, k=2, maxrestarts=20, maxit = 1500)
    plot.mixEM(fit, whichplots = 2, main2=paste0('fit distribution of expression for ', colnames(cpm)[n]))
    mus = fit$mu;
    expressed[jj1, n] = fit$posterior[, which(mus==max(mus))]>posterior.cutoff;
    expressed[jj2, n] = FALSE
  }
  
  par(mfrow=c(1,1))
  EEs = apply(expressed, 1, sum)
  return(EEs>=(ncol(expressed)*fraction.detected.samples))
}

index.outliers = function(data.xx)
{
  c = 1.5
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
  Q1 = quantile(data.xx, 0.25, type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx<lower|data.xx>upper)
  #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
}

Check.3p.5p.arm.switch = function(all)
{
  xx = all[order(all$gene), ]
  ggs = sapply(xx$gene, find.mirName);
  ggs = data.frame(ggs, xx$gene, stringsAsFactors = FALSE)
  yy = as.matrix(xx[, -1]);
  yy[which(is.na(yy))] = 0
  
  Compare.5p.3p = function(x){
    x = as.numeric(x);
    if(x[1]>x[2]) return(1)
    if(x[1]==x[2]) return(0)
    if(x[1]<x[2]) return(-1)
  }
  ggs.uniq = unique(ggs[, 1])
  counts.ps = c()
  names = c()
  for(n in 1:length(ggs.uniq)){
    #n = 1
    jj = which(ggs[, 1]==ggs.uniq[n])
    if(length(jj)>1){
      names = c(names, ggs.uniq[n])
      test = yy[jj, ]
      #means = apply(test, 1, mean)
      ps = apply(test, 2, Compare.5p.3p)
      counts.ps = rbind(counts.ps, c(length(which(ps==1)), length(which(ps==0)), length(which(ps==(-1)))));
    }
  }
  colnames(counts.ps) = c("3p", "3/5p", "5p")
  rownames(counts.ps) = names
  
  length(which(counts.ps[,1]==0 | counts.ps[,3]==0))
  length(which(counts.ps[,1]<=1 | counts.ps[,3]<=1))
  
  mm = match(rownames(counts.ps), list.expressed)
  mm = which(!is.na(mm))
  
  length(which(counts.ps[mm,1]<=1 | counts.ps[mm,3]<=1))
  length(which(counts.ps[mm,1]==0 | counts.ps[mm,3]==0))
}

process.countTable = function(all, design)
{
  index = c()
  for(n in 1:nrow(design))
  {
    #n = 1;
    jj = intersect(grep(design$SampleID[n], colnames(all)), grep("Total.count", colnames(all)));
    if(length(jj)==1) {
      index = c(index,jj)
    }else{print(paste0("ERROR for sample--", design$SampleID[n]))}
  }
  
  newall = data.frame(as.character(all[,1]),  as.matrix(all[, index]), stringsAsFactors = FALSE)
  colnames(newall)[1] = "gene";
  colnames(newall)[-1] = paste0(design$genotype, "_", design$tissue.cell, "_", design$SampleID)
  
  return(newall)
}

cat.countTable = function(xlist)
{
  ## input is list.files for count tables (including directories and file names)
  counts = NULL
  for(n in 1:length(xlist)){
    # n = 1
    ff = read.delim(xlist[n], sep='\t', header = TRUE, as.is = c(1));
    if(n==1){
      ggs = unique(ff[, 1]);
      counts = data.frame(ggs, ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
    }else{
      ggs = unique(c(counts[, 1],ff[, 1]));
      counts = data.frame(ggs, counts[match(ggs, counts[, 1]), -1], ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
    }
  };
  
  colnames(counts)[1] = 'gene'
  return(counts)
  
}

Compare.total.median.normalization = function()
{
  TEST.median.total.normalization = FALSE
  if(TEST.median.total.normalization)
  {
    #cat("nb of expressed miRNAs --", nrow(dds.expressed) )
    cat("size factor is -- ", sizeFactors(dds.expressed), "\n")
    ss = apply(counts(dds.expressed.mature), 2, sum)
    cat("size factor by total expressed mature read counts --- ", ss/mean(ss), "\n")
    ss = apply(counts(dds.expressed), 2, sum)
    cat("size factor by total expressed read counts --- ", ss/mean(ss), "\n")
    #dds <- estimateSizeFactors(dds)
    
    pdfname = paste0(resDir, "Comparision_total_median_normalization_expressed_mature_", specifity, ".pdf") #save all plots during data processing
    pdf(pdfname, width = 12, height = 6)
    cols = rep("blue", nrow(dds.expressed.mature));
    cols[grep("mmu-", rownames(dds.expressed.mature))] = 'red'
    par(mfrow=c(1,2))
    
    for(n in 1:6)
    {
      jj.test = c((2*n-1), 2*n)
      cpm = data.frame(fpm(dds.expressed.mature, robust = FALSE))
      plot(cpm[, jj.test], col=cols, log='xy', main='total normalization'); abline(0, 1, lwd=2.0)
      cpm = data.frame(fpm(dds.expressed.mature, robust = TRUE))
      plot(cpm[, jj.test], col=cols, log='xy', main='median normalization'); abline(0, 1, lwd=2.0)
    }
    dev.off()
  }
  
}

calculate.scaling.factors.using.spikeIns = function(dds, concentrations = c(0.5, 2.5, 5.0, 15, 25, 35, 50, 250),index.spikeIn=c(1:8), read.threshold=1, method="Ratio")
{
  if(method == "Ratio"){
    cpm = fpm(dds, robust = FALSE)
    cpm.spikeIn = cpm[index.spikeIn, ]
    counts.spikeIn = assay(dds)[index.spikeIn, ]
    if(nrow(cpm.spikeIn) != length(concentrations)) stop("wrong number of spike in or concentrations !!!")
    
    norms = rep(NA, ncol(cpm.spikeIn))
    for(n in 1:length(norms))
    {
      #n = 1;
      jj = which(counts.spikeIn[, n]>=read.threshold)
      if(length(jj)>=1)
      {
        norms[n] = 1.0/median((concentrations[jj]) / (cpm.spikeIn[jj, n])); 
        ## NOT use log scale any more 
        #if( ! log.scale ){
        #  fit = lm(concentrations[jj] ~ 0 + cpm.spikeIn[jj,n])
        #  norms[n] = fit$coefficients
        #}else{
        #  norms[n] = exp(median(log(concentrations[jj]) - log(cpm.spikeIn[jj, n])));
        #norms[n]
        #}
      }else{
        cat ("NO spike ins detected or pass the threshod of read number \n")
      }
    }
  }
  
  if(method == "DESeq2"){ # test DESeq2
    require('DESeq2')
    dds.spike= dds[kk, ]
    dds.spike <- dds.spike[ rowSums(counts(dds.spike)) > 5, ]
    dds.spike <- estimateSizeFactors(dds.spike)
    norms = sizeFactors(dds.spike);
  }
  
  if(method == "RUVg"){ # test RUVg method
    library(RUVSeq)
    zfGenes = assay(dds);
    filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
    filtered <- zfGenes[filter,]
    genes <- rownames(filtered)[grep("^cel", rownames(filtered))]
    spikes <- rownames(filtered)[grep("^spike", rownames(filtered))]
    
    x <- as.factor(design$genotype)
    set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))
    set
    
    library(RColorBrewer)
    colors <- brewer.pal(3, "Set2")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set <- betweenLaneNormalization(set, which="upper")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set1 <- RUVg(set, spikes, k=2)
    pData(set1)
    plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set1, col=colors[x], cex=1.2)
  }
  
  return(norms)
  
}


