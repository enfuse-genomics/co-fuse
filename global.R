library(gplots)
#
# Helper functions for appending additional information
# from other databases
#

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#
# Append the new column (colname) to the input data (dat)
# the value in the new column is set to 'geneA','geneB' or'Both'
# if the input gene matches the values in dat$geneA, dat$geneB or both
#
helper.oneGene <- function(dat,colname,gene) {
  gene1 <- trim(gene)
  gene1 <- data.frame(gene=gene1[!duplicated(gene1)],Indicator=1,stringsAsFactors = F)
  dat1 <- left_join(dat,gene1,by=c('geneA'='gene'))
  dat2 <- left_join(dat,gene1,by=c('geneB'='gene'))
  indicator <- cbind(dat1$Indicator,dat2$Indicator)
  indicator[is.na(indicator)] <- 0
  combined <- data.frame(key=2*indicator[,1] + indicator[,2])
  lookupTable <- data.frame(key=0:3,value=c('','GeneB','GeneA','Both'),stringsAsFactors = F)
  combined <- left_join(combined,lookupTable,by=c('key'))
  res <- dat; res[,colname] <- combined$value
  res
}

#
# Append the new column (colname) to the input data (dat)
# the value in the new column is set to TRUE if the input geneA and geneB
# match the values in dat$geneA and dat$geneB
#
helper.twoGenes <- function(dat,colname,geneA,geneB) {
  genes1 <- data.frame(geneA=trim(geneA),geneB=trim(geneB),Indicator=1,stringsAsFactors=F)
  genes1 <- genes1[!duplicated(genes1),]
  dat1 <- left_join(dat,genes1,by=c('geneA'='geneA','geneB'='geneB'))
  dat2 <- left_join(dat,genes1,by=c('geneA'='geneB','geneB'='geneA'))
  indicator <- cbind(dat1$Indicator,dat2$Indicator)
  indicator[is.na(indicator)] <- 0
  indicator <- as.character(rowSums(indicator))
  indicator[indicator=='0'] = ''
  indicator[indicator!=''] = 'Y'
  res <- dat; res[,colname] <- indicator
  res
}

#
# Given the database (DBName), it appends an additional column to 
# indicate if geneA and geneB can be found in the given database
#
lookup.DB <- function(dat,DBName) {
  switch(DBName,
         Human_TSG={
           db <- read.csv('./databases/Human_TSGs.txt',sep='\t',stringsAsFactors = F)
           gene <- db$GeneSymbol
           return(helper.oneGene(dat,'Human_TSG',gene))
         },
         PanCancer_TSG={
           db <- read.csv('./databases/All_down_exp_TSGs_pan-cancer.txt',sep='\t',stringsAsFactors = F)
           gene <- db$GeneName  
           return(helper.oneGene(dat,'PanCaner_TSG',gene))
         },
         Oncogene={
           db <- read.csv('./databases/allOnco_Feb2017.tsv',sep='\t',stringsAsFactors = F)
           gene <- db$symbol
           return(helper.oneGene(dat,'Onco',gene))
         },
         COSMIC={
           db <- read.csv('./databases/Census_allSat Mar 11 09_32_39 2017.csv',stringsAsFactors = F)
           gene <- db$Gene.Symbol
           return(helper.oneGene(dat,'COSMIC',gene))
         },
         Kinase={
           db <- read.csv('./databases/kinase database.csv',stringsAsFactors = F)
           gene <- db$Symbol
           return(helper.oneGene(dat,'Kinase',gene))
         },
         ChimerKB={
           db <- read.csv('./databases/ChimerDB3.0_ChimerKB.csv',stringsAsFactors = F)
           geneA <- db$H_gene; geneB <- db$T_gene
           return(helper.twoGenes(dat,'ChimerKB',geneA,geneB))
         },
         ChimerPub={
           db <- read.csv('./databases/ChimerDB3.0_ChimerPub.csv',stringsAsFactors = F)
           geneA <- db$H_gene; geneB <- db$T_gene
           return(helper.twoGenes(dat,'ChimerPub',geneA,geneB))
         },
         ChimeSeq={
           db <- read.csv('./databases/ChimerDB3.0_ChimerSeq.csv',stringsAsFactors = F)
           geneA <- db$H_gene; geneB <- db$T_gene
           return(helper.twoGenes(dat,'ChimerSeq',geneA,geneB))
         },
         {
           # default, do nothing
           return(dat)
         }
  )
}


fusion.plot.heatmat <- function(data, outfilename,exRowCol)
{
  if(missing(exRowCol)) {
    cexRowValue <- 0.4
    cexColValue <- 0.5
  } else {
    cexRowValue <- exRowCol[1]
    cexColValue <- exRowCol[2]
  }
  data$fusion <- paste(data$geneA, data$geneB, sep=" : ")
  n <- length(colnames(data))-1
  data.1 <- as.matrix(data[, 3:n])
  rownames(data.1) <- data$fusion
  my_palette <- colorRampPalette(c("grey", "black", "red"))
  pdfheight <- ceiling(log(nrow(data.1))*4) + 1
  pdfwidth <- ceiling(log(ncol(data.1))*3) + 1
  pdf(file=outfilename,width=pdfwidth,height=pdfheight,paper='special')
  heatmap.2(data.1, col=my_palette, distfun=function(x) dist(x, method="binary"), 
            hclustfun=function(x) hclust(x, method="ward"), scale="none", 
            cexRow=cexRowValue, cexCol=cexColValue, density.info="none", 
            trace="none", key=TRUE, symkey=FALSE,
            margins=c(8,10))
  dev.off()
}

