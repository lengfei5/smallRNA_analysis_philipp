#### Functions
process.countTable = function(all, design, select.Total.count = TRUE)
{
  newall = data.frame(as.character(all[,1]), stringsAsFactors = FALSE)
  
  for(n in 1:nrow(design))
  {
    #n = 1;
    ## found the matching ID in design matrix
    if(select.Total.count){
      jj = intersect(grep(design$SampleID[n], colnames(all)), grep("Total.count", colnames(all)));
    }else{
      jj = grep(design$SampleID[n], colnames(all));
    }
    
    ## deal with several mapping 
    if(length(jj)==1) {
      #index = c(index,jj)
      newall = data.frame(newall, all[, jj])
    }else{
      cat(length(jj), " samples found for ID", design$SampleID[n], "\n")
      cat("start to merge those samples considering them as technical replicates...\n")
      newall = data.frame(newall, apply(as.matrix(all[, jj]), 1, sum))
    }
  }
  
  colnames(newall)[1] = "gene";
  colnames(newall)[-1] = paste0(design$stage, "_", design$treatment, "_", design$SampleID)
  
  return(newall)
}

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

find.expressed.mature.miRNA.for.dm = function(dds, cpm.threshold=10, sampletoUse="untreated")
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

my.cpm.normalization = function(countData)
{
  cpm = matrix(NA, nrow = nrow(countData), ncol = ncol(countData));
  colnames(cpm) = colnames(countData)
  rownames(cpm) = rownames(countData)
  for(n in 1:ncol(countData))
  {
    cpm[,n] = countData[,n]/sum(countData[,n])*10^6
  }
  
  return(cpm)
}

calculate.scaling.factors.using.spikeIns = function(countData, concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25), index.spikeIn=c(1:8), 
                                                    read.threshold=5, 
                                                    method="Ratio", plot.concentration.vs.cpm = TRUE)
{
  # countData = raw;
  if(method == "Ratio"){
    ss = apply(as.matrix(countData), 2, sum);
    cpm = my.cpm.normalization(countData)
    cpm.spikeIn = cpm[index.spikeIn, ]
    res = matrix(NA, ncol = ncol(cpm), nrow = nrow(cpm))
    rownames(res) = rownames(countData);
    colnames(res) = colnames(countData);
    counts.spikeIn = countData[index.spikeIn, ]
    
    if(nrow(cpm.spikeIn) != length(concentrations)) stop("wrong number of spike in or concentrations !!!")
    
    scaling.factors = rep(NA, ncol(cpm.spikeIn))
    for(n in 1:length(scaling.factors))
    {
      #n = 1;
      jj = which(counts.spikeIn[, n]>=read.threshold)
      if(length(jj)>=1)
      {
        scaling.factors[n] = median((concentrations[jj]) / (cpm.spikeIn[jj, n])); 
        res[,n] = cpm[,n]*scaling.factors[n];
        
        if(plot.concentration.vs.cpm){
          #kk = grep("spikeIn", rownames(cpm))
          plot((cpm.spikeIn[, n]+10^-6),  concentrations, log='xy', 
               xlab="cpm", ylab="amol per mug total RNA", main=paste0(colnames(cpm)[n], '--', signif(scaling.factors[n], d=2)), cex=3.0, col='darkgreen', pch= 16);
          abline(h=10^-6, lwd=2.0, col='darkgray')
          points(range((cpm.spikeIn[, n]+10^-6)), range((cpm.spikeIn[, n]+10^-6))*scaling.factors[n], type = 'l', lwd=3.0, col='darkblue', lty=1)
        }
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
    #library(zebrafishRNASeq)
    #data("zfGenes")
    #countData= raw;
    #zfGenes = assay(dds);
    #filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
    #filtered <- zfGenes[filter,]
    genes <- rownames(countData)[grep("^cel", rownames(countData))]
    spikes <- rownames(countData)[grep("^spike", rownames(countData))]
    
    x <- as.factor(design$genotype)
    set <- newSeqExpressionSet(as.matrix(countData), phenoData = data.frame(x, row.names=colnames(countData)))
    set
    
    library(RColorBrewer)
    par(mfrow=c(3,2))
    colors <- brewer.pal(3, "Set2")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set <- betweenLaneNormalization(set, which="upper")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set1 <- RUVg(set, spikes, k=1)
    pData(set1)
    plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set1, col=colors[x], cex=1.2)
    
    cpm.ruv = my.cpm.normalization(normCounts(set1))
  }
  
  ## compute the size factor for DESeq2
  ss.geo.mean = prod(ss)^(1/length(ss))
  norms = 1/scaling.factors * ss / ss.geo.mean
  
  return(list(scaling.factors = scaling.factors, norms4DESeq2 = norms, cpm = cpm, normalization.spikeIn = res))
  
}

Merge.techinical.replicates.for.N2 = function(all, id.list=list(c("57751", "57753"), c("57752", "57754")))
{
  for(n in 1:length(rep.technical))
  {
    index = c()
    for(id in rep.technical[[n]])
    {
      #print(id)
      index = c(index, which(design$SampleID==id))
    }
    
    design$SampleID[index[1]] = paste0(design$SampleID[index], collapse = ".")
    ss = apply(all[, (index+1)], 1, function(x) sum(x, na.rm = TRUE))
    all[, (index[1]+1)] = ss;
    colnames(all)[(index[1]+1)] = paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  design$SampleID[index[1]])
    design = design[-index[-1], ]
    all = all[, -(index[-1]+1)]
  }
}

identify.expressed.miRNAs= function(countData, design.matrix, cpm.threshold=10) 
{
  cat('identify list of expressed miRNAs for each stage\n')
  require('DESeq2')
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ treatment + stage)
  cpm =  fpm(dds, robust = FALSE)
  
  colnames(cpm) = paste0(colnames(cpm), ".cpm")
  
  test.theshold = FALSE
  if(test.theshold){
    par(mfrow=c(4,4))
    for(n in 1:ncol(cpm))
    {
      hist(log2(cpm[which(cpm[,n]>0), n]), breaks = 30,  xlab="log2(cpm)", main=colnames(cpm)[n], cex=3.0, col='darkgray', pch= 16);
      ggs = rownames(cpm)[which(cpm[,n]>10)]
      ggs = ggs[grep('spikeIn_', ggs, invert = TRUE)]
      ggs = sapply(ggs, function(x) gsub("-3p", '', x))
      ggs = sapply(ggs, function(x) gsub("-5p", '', x))
      cat("nb of miRNAs above the threshold--", length(unique(ggs)), "\n")
      abline(v=log2(10), lwd=2.0, col='red')
      #points(range((cpm[kk, n]+10^-6)), range((cpm[kk, n]+10^-6))/norms[n], type = 'l', lwd=3.0, col='darkblue', lty=1)
    }
  }
  
  cpm.untreated = cpm[, grep("_untreated", colnames(cpm))]
  kk.spikeIn = grep("spikeIn_", rownames(cpm.untreated))
  if(length(kk.spikeIn)>0) cpm.untreated = cpm.untreated[-kk.spikeIn, ]
  
  times = unique(design.matrix$stage[which(design.matrix$treatment=="untreated")])
  compare.threshold = c()
  for(n in 1:nrow(cpm.untreated))
  {
    test = 0;
    for(t in times) 
    {
      #cat(test, "---", mean(cpm.untreated[n, grep(t, colnames(cpm.untreated))]))
      if(mean(cpm.untreated[n, grep(t, colnames(cpm.untreated))])>cpm.threshold)  test = test + 1;
    }
    compare.threshold = c(compare.threshold, test)
  }
  
  expressed = data.frame(rownames(cpm.untreated), sapply(rownames(cpm.untreated), find.mirName), compare.threshold, cpm.untreated, stringsAsFactors = FALSE)
  colnames(expressed)[1:3] = c('miRNA', 'gene', 'nb.stage.above.threshold' )
  gg.expressed = unique(expressed$gene[which(expressed$nb.stage.above.threshold>0)])
  
  mm = match(expressed$gene, gg.expressed)
  expressed = expressed[which(!is.na(mm)==TRUE), ]
  
  return(expressed)
  
}

find.mature.ones.for.expressed.miRNAs = function(list.expressed.miRNAs)
{
  list.expressed = list.expressed.miRNAs
  colnames(list.expressed)[c(1,2)] = c("miRNA", "gene")
  list.expressed = list.expressed[order(list.expressed$miRNA), ]
  ggs.uniq = unique(list.expressed$gene)
  
  mature = rep(NA, nrow(list.expressed))
  kk = grep('.cpm', colnames(list.expressed))
  
  for(n in 1:length(ggs.uniq))
  {
    jj = which(list.expressed$gene==ggs.uniq[n])
    if(length(jj)>1){
      index.max = apply(list.expressed[jj, kk], 2, which.max)
      nb.first.max = length(which(index.max==1))
      nb.sec.max = length(which(index.max==2))
      #cat(n, ": ", as.character(ggs.uniq[n]), "--",  nb.first.max, "--", nb.sec.max, "\n")
      if(nb.first.max>nb.sec.max){
        mature[jj[1]] = TRUE; mature[jj[2]] = FALSE;  
      }else{
        mature[jj[1]] = FALSE; mature[jj[2]] = TRUE; 
      }
    }else{
      #cat(n,": ", as.character(ggs.uniq[n]),  "-- no selection \t")
      mature[jj] = TRUE;
    }
  }
  expressed.miRNAs = data.frame((list.expressed[, c(1, 2)]), mature=(mature), list.expressed[, -c(1:2)], stringsAsFactors = FALSE)
  return(expressed.miRNAs)
}

find.replicates.by.removing.ID = function(x)
{
  infos = unlist(strsplit(x, "_"))
  return(paste0(infos[-length(infos)], collapse = '_'))
}

average.biological.replicates = function(cpm)
{
  # cpm = cpm.piRNA.bc;
  samples = sapply(colnames(cpm), find.replicates.by.removing.ID, USE.NAMES = FALSE)
  samples.uniq = unique(samples)
  
  if(length(samples.uniq) == length(samples)){
    cat('---no replicates exist---')
  }else{
    cpm.mean = matrix(NA, nrow = nrow(cpm), ncol=length(samples.uniq))
    rownames(cpm.mean) = rownames(cpm)
    colnames(cpm.mean) = samples.uniq;
    
    for(n in 1:ncol(cpm.mean))
    {
      kk = which(samples == colnames(cpm.mean)[n])
      if(length(kk)>1){
        cpm.mean[ ,n] = apply(as.matrix(cpm[,kk]), 1, mean)
      }else{
        if(length(kk)==1) cpm.mean[, n] = cpm[,kk]
      }
    }
  }
  
  return(cpm.mean)
  
}

remove.batch.using.N2.untreated = function(cpm, design.matrix, method = "linear.model")
{
  # cpm = cpm.piRNA
  logcpm = log2(cpm + 2^-6)
  
  if(method == 'linear.model'){
    ## use the linear model for data the log2 scale (sample by sample, NOT gene by gene) to remove the batch effect; use N2 untreated condtion as references
    ## (NOT used here)
    cat('Warnings -- use limma or ComBat instead of this method ')
    
    #reference = "N2_whole.body_untreated"
    #genotypes = sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[1], USE.NAMES = FALSE)
    #ns =  sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
    #treatment = sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[3], USE.NAMES = FALSE)
    pheno = data.frame(design.matrix, stringsAsFactors = FALSE)
    
    jj.N2.untreated = which(pheno$genotype=="N2" & pheno$treatment == "untreated" & pheno$batch == 1)
    batchs = unique(pheno$batch)
    batchs = batchs[which(batchs != 1)]
    for(n in 1:length(batchs))
    {
      jj.untreated = which(pheno$batch == batchs[n] & pheno$treatment == "untreated")
      if(length(jj.untreated)>0)
      {
        if(length(jj.untreated) == 1) x = logcpm[, jj.untreated]
        if(length(jj.untreated)>1) x = apply(logcpm[, jj.untreated], 1, mean)
        y = apply(logcpm[, jj.N2.untreated], 1, mean)
        fit = lm(y ~ x)
        #plot(x, y);
        #abline(fit, col='red', lwd=2.0)
        jj.to.correct = which(pheno$batch == batchs[n])
        for(j in jj.to.correct) {
          logcpm[, j]  = logcpm[, j]*fit$coefficients[2] + fit$coefficients[1];
          cat(colnames(cpm)[j], "-intercept-", fit$coefficients[1], " - slop -", fit$coefficient[2],  "\n")   
        }
      }
    }
    logcpm.bc = logcpm;
  }
  
  if(method == 'limma'){
    
    cat('remove the batch effect using limma \n')
    require('limma')
    design.tokeep = design.matrix
    design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    design.tokeep<-model.matrix(~0 + tissue.cell,  data = design.tokeep)
    logcpm.bc = removeBatchEffect(logcpm, batch = design.matrix$batch, design = design.tokeep)
    
  }
  if(method == 'combat'){
    ## here we use the combat to remove the batch effect 
    ## the combat requires the log2cpm
    cat('remove the batch effect using ComBat \n')
    require("sva")
    # example from the ComBat function in the R package 'sva'
    TEST.example = FALSE
    if(TEST.example){
      library(bladderbatch)
      data(bladderdata)
      dat <- bladderEset[1:50,]
      pheno = pData(dat)
      edata = exprs(dat)
      batch = pheno$batch
      mod = model.matrix(~as.factor(cancer), data=pheno)
      combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
    }
    
    batch = design.matrix$batch;
    design.tokeep = design.matrix
    design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    mod = model.matrix(~ as.factor(tissue.cell), data = design.tokeep);
    #conds = data.frame(rep(c("untreated", "treated"), ncol(cpm)/2))
    #colnames(conds) = 'treatment'
    #mod = model.matrix(~ as.factor(treatment), conds)
    logcpm.bc = ComBat(dat=logcpm, batch=batch, mod=mod, par.prior=TRUE, ref.batch=1)
  }
  
  return(2^logcpm.bc)
}

Test.piRNA.normalization.batch.removal = function(cpm, design.matrix)
{
  ## here is a function to test piRNA normalization and batchRemoval for the gene expression matrix in deconvolution analysis
  ## there will be one PCA plot and four samll plots for the positive controls
  # cpm = cpm.piRNA
  main.names = deparse(substitute(cpm)) 
  if(any(cpm==0)) cpm = cpm + 2^-6 
  
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr");
  
  ## pca plots only for untreated samples
  pca = prcomp(t(log2(cpm)), scale. = FALSE, center = FALSE)
  pca2save = data.frame(pca$x, condition=design.matrix$treatment, batch = design.matrix$batch, name=colnames(cpm), tissue = design.matrix$tissue.cell)
  sels = which(pca2save$condition == "untreated")
  ggp = ggplot(data=pca2save[sels, ], aes(PC1, PC2, label = batch, color = tissue)) + geom_point(size=4) +
    geom_text(hjust = 0.1, nudge_y = 0.1, size=5) +
    ggtitle(paste0("PCA - ", main.names))
  plot(ggp);
  #pca=plotPCA(cpm, intgroup = colnames(design.matrix)[c(3, 5, 7)], returnData = FALSE)
  #print(pca)
  
  ## ratio between treated and untreated to test if piRNA normalization makes sense
  ns = unique(design.matrix$tissue.cell)
  ratios = matrix(NA, ncol = length(ns), nrow = nrow(cpm))
  rownames(ratios) = rownames(cpm);
  colnames(ratios) = ns
  total = ratios;
  tcs = ratios;
  for(n in 1:length(ns))
  {
    #jj = which(design.matrix$tissue.cell==ns[n])
    total[,n] = apply(cpm[, which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="untreated")], 1, mean)
    tcs[,n] = apply(cpm[, which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="treated")], 1, mean)
    ratios[,n] = tcs[,n]/total[,n] 
  }
  
  par(mfrow=c(2, 2))
  hist(log10(ratios), xlab = "log10(treated/untreated)", main = paste0(main.names), breaks = 100)
  abline(v=c(-3:0), col='darkred', lwd=2.0)
  
  ## check lsy-6 in ASE, Glutatamatergic and Ciliated neurons
  kk = which(rownames(cpm)=='lsy-6')
  lims = range(c(total[kk, ], tcs[kk,]))
  plot(c(1:length(ns)), tcs[kk, ], type= 'b', col='darkblue', cex=2.0, log='y', ylim =lims, main = paste0("lsy-6 in ", main.names), xlab=NA, 
       ylab = 'normalizaed by piRNAs (in log)', axe = FALSE)
  points(c(1:length(ns)), total[kk, ], type= 'b', col='black', cex=2.0)
  legend("topright", col=c('darkblue', "black"),  bty = "n", legend = c("treated", "untreated"), lty=1 )
  axis(2, las= 1)
  ns.short = sapply(ns, function(x) unlist(strsplit(x, "[.]"))[1], USE.NAMES = FALSE)
  ns.short[c(3:5, 8,10, 11)] = c("Seroton", "Dopamin", "Gluta", "mecha", 'pharyn', "cholin")
  axis(1, at=c(1:length(ns)), labels = ns.short, las=2,cex=0.5)
  box()
  
  plot(tcs[,which(colnames(tcs) == "Dopaminergic.neurons")], tcs[,which(colnames(tcs) == "Ciliated.sensory.neurons")], log='xy', 
       xlab='Dopaminergic (in log)', ylab='Ciliated (in log)', main = main.names)
  abline(0, 1, lwd=2.0, col='red')
  
  plot(tcs[,which(colnames(tcs) == "mechanosensory.neurons" )], tcs[,which(colnames(tcs) == "unc-86.expressing.neurons")], log='xy', 
       xlab='mechanosensory (in log)', ylab='unc-86 (in log)', main = main.names)
  abline(0, 1, lwd=2.0, col='red')
  
}

calculate.pvalues.two.groups.overlapping = function(nb.total, nb.group.A, nb.group.B, nb.overlapping)
{
  total = as.numeric(nb.total)
  q = as.numeric(nb.overlapping)
  m = as.numeric(nb.group.A)
  nn = total - m;
  k = as.numeric(nb.group.B)
  rr = q/(k*m/total)
  pvals = phyper((q-1), m, nn, k, lower.tail = FALSE, log.p = FALSE)
  
  return(c(rr, pvals))
}


