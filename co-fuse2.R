#!/usr/bin/env Rscript

#
# Assumptions:
# 1. The first two columns of the software output must be "GeneA" and "GeneB"
#

library(dplyr)

#
# Return fusion with the number of samples 
#
countSamples <- function(inputDir) { 
  folders <- list.dirs(path=inputDir,full.names=TRUE,recursive=TRUE)
  sampleNames <- list.dirs(path=inputDir,full.names=FALSE,recursive=TRUE)
  numFolders <- length(folders)
  
  listAB <- vector(mode="list", length=numFolders)
  
  # for each files in the folder
  for (i in 2:numFolders) {
    # load to data frame
    filename <- list.files(path=folders[i],pattern="*.GRCh37.txt",full.name=TRUE,recursive=TRUE)
    dat <- read.csv(file=filename,sep="\t",stringsAsFactors = F)
    df <- data.frame(geneA=dat[,1],geneB=dat[,2],stringsAsFactors = F)
    # remove duplicates within the file
    listAB[[i]] <- df %>% group_by(geneA,geneB) %>% filter(row_number() == 1)
  }
  
  # union all genes
  merged.AB <- Reduce(function(...) merge(...,all=TRUE),listAB[2:numFolders])
  
  for (i in 2:numFolders) {
    tmpAB <- listAB[[i]]
    if(nrow(tmpAB)>0) { 
      tmpAB[[sampleNames[i]]] <- 1
      merged.AB <- merged.AB %>% left_join(.,tmpAB,by=c("geneA","geneB"))
    } else {
      merged.AB[[sampleNames[i]]] <- 0
    }
  }
  
  # exclude a pair of genes (geneA, geneB) that occur only in one sample
  merged.AB[is.na(merged.AB)] <- 0
  
  #############################################################################
  # Summarize the results (gene, #samples, total samples)
  ############################################################################# 
  res.AB <- data.frame(geneA=merged.AB$geneA,geneB=merged.AB$geneB,
                       NumSamples=0,TotalSamples=numFolders-1,stringsAsFactors = F)
  res.AB$NumSamples <- rowSums(merged.AB[,-(1:2)])
  
  return(res.AB)
}


FisherTest <- function(inputDirGroup1,inputDirGroup2,outputDir) { 
  folders1 <- list.dirs(path=inputDirGroup1,full.names=TRUE,recursive=TRUE)
  folders2 <- list.dirs(path=inputDirGroup2,full.names=TRUE,recursive=TRUE)

  if (length(folders1) < 1) {
    stop("The input folder (first group) must contain at least one sample!!")
  } else if (length(folders2) < 1) {
    stop("The input folder (second group) must contain at least one sample!!")
  }
  

  group1 <- countSamples(inputDirGroup1)
  Num.group1 <- paste0("Num.",basename(inputDirGroup1))
  Total.group1 <- paste0("Total.",basename(inputDirGroup1))
  colnames(group1) <- c("geneA","geneB",Num.group1,Total.group1)

  group2 <- countSamples(inputDirGroup2)
  Num.group2 <- paste0("Num.",basename(inputDirGroup2))
  Total.group2 <- paste0("Total.",basename(inputDirGroup2))
  colnames(group2) <- c("geneA","geneB",Num.group2,Total.group2)
  
  combined <- full_join(group1,group2,by=c("geneA","geneB"))
  combined[is.na(combined)] <- 0
  combined[,Total.group1] <- group1[1,4]
  combined[,Total.group2] <- group2[1,4]
  
  #
  # Filter out non-recurrent gene
  #
  combined <- combined %>% filter_(paste('! (',Num.group1,'==',1,' & ',Num.group2, '==', 0,')') ) %>%
    filter_(paste('! (',Num.group1, '==', 0, '& ',Num.group2, '==', 1,')') ) %>%
    filter_(paste('! (',Num.group1, '==', 1, '& ',Num.group2, '==', 1,')') )  
  

  group1_YN <- sapply(combined[,Num.group1], function(x) ifelse(x>0,"Y","N"))
  group2_YN <- sapply(combined[,Num.group2], function(x) ifelse(x>0,"Y","N"))
  

  tmp <- combined[,c(Num.group1,Total.group1,Num.group2,Total.group2)]
  colnames(tmp) <- c('Num1','Total1','Num2','Total2')
  #
  #                        Group_1              Group_2
  # Fusion_A             Num_Group_1          Num_Group_2
  # Non Fusion_A      Total-Num_Group_1    Total-Num_Group_2
  #
  fisher <- tmp %>% rowwise() %>%
                    mutate(pvalue=fisher.test(matrix(c(Num1,Total1-Num1,
                                                       Num2,Total2-Num2),2,2))$p.value)
  
  combined$pvalue <- fisher$pvalue
  combined[[basename(inputDirGroup1)]] <- group1_YN
  combined[[basename(inputDirGroup2)]] <- group2_YN
  
  res <- combined[,c(1,2,8,9,7,3,4,5,6)]
  
  # sort them based on p-value
  res <- res[with(res,order(pvalue)),]
  
  
  #===========================================================================
  # save results
  #===========================================================================
  # Create output directory if not exist
  dir.create(outputDir,showWarnings = F)
  #############################################################################
  # save the summary results (gene, #samples, total samples, sample_ids)
  #############################################################################
  write.csv(res,file=paste0(outputDir,'/twoGroups.csv'),row.names=F,quote=F)
}


options(scipen=99)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop(paste("Three arguments must be provided.\n",
             "Usage: Rscript co-fuse2.R _GROUP_1_INPUT_FOLDER_  _GROUP_2_INPUT_FOLDER_  _OUTPUT_FOLDER_  \n\n")
             ,call.=FALSE)
} else {
  cat("Evaluating Fisher's Exact test\n")
  FisherTest(args[1],args[2],args[3])
}



