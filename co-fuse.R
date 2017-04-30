#!/usr/bin/env Rscript


#ONCOFUSE_JAR_FILE <- '/windows/D/Fusion/oncofuse-1.1.1/Oncofuse.jar'


#
# Assumptions:
# 1. The first two columns of the software output must be "GeneA" and "GeneB"
#

library(dplyr)
library(tsne)
library(gplots)

source('./global.R')

fusion.plot.heatmat <- function(data, outfilename)
{
  data$fusion <- paste(data$geneA, data$geneB, sep="-")
  n <- length(colnames(data))-1
  data.1 <- as.matrix(data[, 3:n])
  rownames(data.1) <- data$fusion
  my_palette <- colorRampPalette(c("grey", "black", "red"))
  pdfheight <- ceiling(log(nrow(data.1))*4) + 1
  pdfwidth <- ceiling(log(ncol(data.1))*3) + 1
  pdf(file=outfilename,width=pdfwidth,height=pdfheight,paper='special')
  heatmap.2(data.1, col=my_palette, distfun=function(x) dist(x, method="binary"), 
            hclustfun=function(x) hclust(x, method="ward"), scale="none", 
            cexRow=0.4, cexCol=0.5, density.info="none", trace="none", key=TRUE, symkey=FALSE)
  dev.off()
}


fusion.plot.tsne <- function(data, perplexity=perplexity, outfilename)
{
  #data$fusion <- paste(data$geneA, data$geneB, sep="-")
  n <- length(colnames(data))
  data.1 <- as.matrix(data[, 3:n])
  pdfwidth <- ceiling(log(ncol(data.1))*3) + 1
  pdf(file=outfilename,width=pdfwidth,height=pdfwidth,paper='special')
  set.seed(12345)
  d <- dist(t(data.1), method="binary")
  tsne.data <- tsne(d,perplexity=perplexity)
  rownames(tsne.data)<- colnames(data.1)
  plot(tsne.data[,1], tsne.data[,2], type="n", xlab="tsne-x", ylab="tsne-y")
  text(tsne.data[,1], tsne.data[,2],labels=row.names(tsne.data), cex=0.8)
  abline(h=0,v=0,col="gray75")
  dev.off()
}


#
# ======================================================
#





countGene <- function(software,inputDir,outputDir,tsne_perplexity=5) { 
  fusionSoftware <- tolower(software)
  if (fusionSoftware == 'fusioncatcher') {
    filepattern <- '*.GRCh37.txt'
    geneACol <- 1; geneBCol <- 2; headerRow = TRUE;
  } else if (fusionSoftware == 'defuse') {
    filepattern <- '*.filtered.tsv'
    geneACol <- 31; geneBCol <- 32; headerRow = TRUE;
    
  } else if (fusionSoftware == 'tophat') {
    filepattern <- '*.txt'
    geneACol <- 2; geneBCol <- 5; headerRow = FALSE;
  } else if (fusionSoftware == 'soapfuse') {
    filepattern <- '*.specific.for.genes'
    geneACol <- 1; geneBCol <- 6; headerRow = TRUE;
  } else if (fusionSoftware == 'generic') {
    filepattern <- '*.txt'
    geneACol <- 1; geneBCol <- 2; headerRow = TRUE;
  } else {
    cat("For supporting of other softwares, please contact us.\n")
    return
  }
  
  
  folders <- list.dirs(path=inputDir,full.names=TRUE,recursive=TRUE)
  sampleNames <- list.dirs(path=inputDir,full.names=FALSE,recursive=TRUE)
  numFolders <- length(folders)
  
  listAB <- vector(mode="list", length=numFolders)
  listA <- vector(mode="list", length=numFolders)
  listB <- vector(mode="list", length=numFolders)
  
  softwareOutputs <- data.frame()
  
  # for each files in the folder
  for (i in 2:numFolders) {
    # load to data frame
    filename <- list.files(path=folders[i],pattern=filepattern,full.name=TRUE,recursive=TRUE)
    dat <- read.csv(file=filename,sep="\t",header=headerRow,stringsAsFactors = F)
    df <- data.frame(geneA=dat[,geneACol],geneB=dat[,geneBCol],stringsAsFactors = F)
    # remove duplicates within the file
    listAB[[i]] <- df %>% group_by(geneA,geneB) %>% filter(row_number() == 1)
    listA[[i]]  <- df %>% select(geneA) %>% group_by(geneA) %>% filter(row_number() == 1)
    listB[[i]]  <- df %>% select(geneB) %>% group_by(geneB) %>% filter(row_number() == 1)
    # save original software output
    if (nrow(dat) > 0) {
      softwareOutputs <- rbind(softwareOutputs,data.frame(dat,Sample=sampleNames[i]))
    }
  }
  
  # union all genes
  merged.AB <- Reduce(function(...) merge(...,all=TRUE),listAB[2:numFolders])
  merged.A <- Reduce(function(...) merge(...,all=TRUE),listA[2:numFolders])
  merged.B <- Reduce(function(...) merge(...,all=TRUE),listB[2:numFolders])
  
  for (i in 2:numFolders) {
    tmpAB <- listAB[[i]]
    if(nrow(tmpAB)>0) { 
      tmpAB[[sampleNames[i]]] <- 1
      merged.AB <- merged.AB %>% left_join(.,tmpAB,by=c("geneA","geneB"))
    } else {
      merged.AB[[sampleNames[i]]] <- 0
    }
    tmpA <- listA[[i]]
    if(nrow(tmpA)>0) { 
      tmpA[[sampleNames[i]]] <- 1 
      merged.A <- merged.A %>% left_join(.,tmpA,by=c("geneA"))
    } else {
      merged.A[[sampleNames[i]]] <- 0
    }
    tmpB <- listB[[i]]; 
    if(nrow(tmpB)>0) { 
      tmpB[[sampleNames[i]]] <- 1 
      merged.B <- merged.B %>% left_join(.,tmpB,by=c("geneB"))
    } else {
      merged.B[[sampleNames[i]]] <- 0
    }
  }
  
  # exclude a pair of genes (geneA, geneB) that occur only in one sample
  merged.AB[is.na(merged.AB)] <- 0
  numSamples.AB <- rowSums(merged.AB[,3:ncol(merged.AB)])
  merged.AB <- merged.AB[numSamples.AB>1,]
  
  merged.A[is.na(merged.A)] <- 0
  numSamples.A <- rowSums(merged.A[,2:ncol(merged.A)])
  merged.A <- merged.A[numSamples.A>1,]
  
  merged.B[is.na(merged.B)] <- 0
  numSamples.B <- rowSums(merged.B[,2:ncol(merged.B)])
  merged.B <- merged.B[numSamples.B>1,]
    
  
  #############################################################################
  # Summarize the results (gene, #samples, total samples, sample_ids)
  ############################################################################# 
  res.AB <- data.frame(geneA=merged.AB$geneA,geneB=merged.AB$geneB,
                       NumSamples=0,TotalSamples=numFolders-1,SampleIDs="",stringsAsFactors = F)
  res.AB$NumSamples <- rowSums(merged.AB[,-(1:2)])
  for (i in 1:nrow(merged.AB)) { 
    res.AB$SampleIDs[i] <- paste(colnames(merged.AB)[colSums(merged.AB[i,]==1)>0],collapse=';')
  }
  res.AB <- res.AB[with(res.AB,order(-NumSamples)),] # sort by Number of Samples
  
  res.A <- data.frame(geneA=merged.A$geneA,
                      NumSamples=0,TotalSamples=numFolders-1,SampleIDs="",stringsAsFactors = F)
  res.A$NumSamples <- rowSums(merged.A[,-1])
  for (i in 1:nrow(merged.A)) {
    res.A$SampleIDs[i] <- paste(colnames(merged.A)[colSums(merged.A[i,]==1)>0],collapse=';')
  }
  res.A <- res.A[with(res.A,order(-NumSamples)),] # Sort by Number of Samples
  
  res.B <- data.frame(geneB=merged.B$geneB,
                      NumSamples=0,TotalSamples=numFolders-1,SampleIDs="",stringsAsFactors = F)
  res.B$NumSamples <- rowSums(merged.B[,-1])
  for (i in 1:nrow(merged.B)) {
    res.B$SampleIDs[i] <- paste(colnames(merged.B)[colSums(merged.B[i,]==1)>0],collapse=';')
  }
  res.B <- res.B[with(res.B,order(-NumSamples)),] # Sort by Number of Samples
  
  
  # For each geneA, list geneB that can pair up with this geneA
  # also count the number of geneBs
  res.A1 <- res.AB %>% mutate(geneB=paste0(geneB,'(',NumSamples,')')) %>%
            select(geneA,geneB) %>% group_by(geneA) %>%
            summarize(numGeneBs=n(),geneBs=paste(geneB, collapse=';'))
  res.A1 <- res.A1[with(res.A1,order(-numGeneBs)),]
  
  res.B1 <- res.AB %>% mutate(geneA=paste0(geneA,'(',NumSamples,')')) %>%
            select(geneA,geneB) %>% group_by(geneB) %>% 
            summarize(numGeneAs=n(),geneAs=paste(geneA, collapse=';'))
  res.B1 <- res.B1[with(res.B1,order(-numGeneAs)),]
  
  
  
  
  #===========================================================================
  # save results
  #===========================================================================
  # Create output directory if not exist
  dir.create(outputDir,showWarnings = F)
  #############################################################################
  # save the summary results (gene, #samples, total samples, sample_ids)
  #############################################################################
  write.csv(res.AB,file=paste0(outputDir,'/summary.AB.csv'),row.names=F,quote=F)
  write.csv(res.A,file=paste0(outputDir,'/summary.A.csv'),row.names=F,quote=F)
  write.csv(res.B,file=paste0(outputDir,'/summary.B.csv'),row.names=F,quote=F)
  
  write.csv(res.A1,file=paste0(outputDir,'/PossiblePair.A.csv'),row.names=F,quote=F)
  write.csv(res.B1,file=paste0(outputDir,'/PossiblePair.B.csv'),row.names=F,quote=F)
  
  
  #############################################################################
  # save original software output into one file
  #############################################################################
  #tmp.softwareOutputs <- softwareOutputs
  #tmp.res.AB <- res.AB %>% select(geneA,geneB)
  #colnames(tmp.softwareOutputs)[c(geneACol,geneBCol)] <- c("geneA","geneB")
  #res.keepFusion <- tmp.softwareOutputs %>% inner_join(tmp.res.AB, by=c("geneA","geneB"))
  #colnames(res.keepFusion)[c(geneACol,geneBCol)] <- colnames(softwareOutputs)[c(geneACol,geneBCol)]
  #write.table(res.keepFusion,file=paste0(outputDir,'/OriginalSoftwareOutput.txt'),row.names=F,quote=F,sep='\t')
  
  
  tmp.softwareOutputs <- softwareOutputs
  tmp.res.AB <- res.AB %>% select(geneA,geneB,NumSamples)
  
  # lookup geneA and geneB in other databases
  DB <- c('Human_TSG', 'PanCancer_TSG', 'Oncogene', 'COSMIC', 'Kinase',
          'ChimerKB', 'ChimerPub', 'ChimeSeq'
          )
  for (db in DB) {
    tmp.res.AB <- lookup.DB(tmp.res.AB,db)
  }
  
  
  colnames(tmp.softwareOutputs)[c(geneACol,geneBCol)] <- c("geneA","geneB")
  res.keepFusion <- inner_join(tmp.softwareOutputs,tmp.res.AB, by=c("geneA","geneB"))
  # remove duplicate
  #res.keepFusion <- res.keepFusion %>% group_by(geneA,geneB,
  #                                              Fusion_point_for_gene_1.5end_fusion_partner.,
  #                                              Fusion_point_for_gene_2.3end_fusion_partner.
  #) %>% filter(row_number() == 1)
  res.keepFusion <- res.keepFusion[with(res.keepFusion,order(-NumSamples)),] # Sort by Number of Samples
  res.keepFusion$NumSamples <- NULL
  colnames(res.keepFusion)[c(geneACol,geneBCol)] <- colnames(softwareOutputs)[c(geneACol,geneBCol)]
  write.table(res.keepFusion,file=paste0(outputDir,'/OriginalSoftwareOutput.txt'),row.names=F,quote=F,sep='\t')
  
  
  # Run Oncofuse
  #tissue_types <- c('HEM','EPI','MES','AVG')
  #for (t in tissue_types) {
  #  cmd <- paste0('java -Xmx2G -jar ',ONCOFUSE_JAR_FILE,' ',
  #                outputDir,'/OriginalSoftwareOutput.txt ',
  #                'fcatcher ',t,' ',
  #                outputDir,'/oncofuse_output_',t,'.txt')
  #  print(cmd)
  #  system(cmd)
  #}
  
  
  
  
  #############################################################################
  # save figure (bar plot)
  #############################################################################  
  barplot.AB.filename <- paste0(outputDir,'/barplot.AB.pdf')
  pdfheight <- ceiling(log(nrow(res.AB))*4) + 1
  pdf(barplot.AB.filename,width=7,height=pdfheight,paper='special')
  par(las=2)
  par(mar=c(5,8,4,2)) # bottom,left,top,right
  yaxisName <- res.AB %>% mutate(geneAB=paste0(geneA,';',geneB)) %>% select(geneAB)
  barplot(res.AB$NumSamples/res.AB$TotalSamples,main="Sample Distribution",horiz=TRUE,
          names.arg=yaxisName$geneAB, xlab="Sample Ratio", cex.axis=0.6,cex.names = 0.4)
  dev.off()
  
  
  
  #############################################################################
  # save tables (1 and 0)
  #############################################################################
  #write.csv(merged.AB,file=paste0(outputDir,'/table.AB.csv'),row.names=F,quote=F)
  #write.csv(merged.A,file=paste0(outputDir,'/table.A.csv'),row.names=F,quote=F)
  #write.csv(merged.B,file=paste0(outputDir,'/table.B.csv'),row.names=F,quote=F)
  
  #############################################################################
  # Plot the heat map
  #############################################################################  
  heatmap.AB.filename <- paste0(outputDir,'/heatmap.pdf')
  fusion.plot.heatmat(merged.AB, heatmap.AB.filename)
  
  #############################################################################
  # Plot tsne
  #############################################################################   
  tsne.AB.filename <- paste0(outputDir,'/tsne_plot.pdf')
  fusion.plot.tsne(merged.AB, perplexity=tsne_perplexity, tsne.AB.filename)
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop(paste("\n\n--------------------------------------------------------------------------------------\n\n",
             "Four arguments must be provided.\n",
             "Usage: Rscript co-fuse.R _SOFTWARE_ _INPUT_FOLDER_  _OUTPUT_FOLDER_  _TSNE_PERPLEXITY_\n\n",
             "where _SOFTWARE_ is one of 'defuse', 'fusioncatcher', 'tophat', 'soapfuse' or 'generic'\n",
             "      _INPUT_FOLDER_ is the folder containing the output of fusion software\n",
             "      _OUTPUT_FOLDER_ is the folder containing the output of Co-Fuse\n",
             "      _TSNE_PERPLEXITY_ is the parameter for TSNE (default 5)\n",
             "\n\n--------------------------------------------------------------------------------------\n\n")
             ,call.=FALSE)
} else {
  cat("Evaluating the data set\n")
  countGene(args[1],args[2],args[3],as.numeric(args[4]))
}


