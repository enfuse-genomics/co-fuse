#!/usr/bin/env Rscript

#
# Assumptions:
# 1. The first two columns of the software output must be "GeneA" and "GeneB"
#

library(dplyr)
library(rpart)
library(rpart.plot)

source('./global.R')

#
# Extract a table where each row represents geneA/geneB
# and each col represents samples
#
extractFeatures <- function(software,inputDir) { 
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
  
  
  # for each files in the folder
  for (i in 2:numFolders) {
    # load to data frame
    filename <- list.files(path=folders[i],pattern=filepattern,full.name=TRUE,recursive=TRUE)
    dat <- read.csv(file=filename,sep="\t",header=headerRow,stringsAsFactors = F)
    df <- data.frame(geneA=dat[,geneACol],geneB=dat[,geneBCol],stringsAsFactors = F)
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
  
  return(merged.AB)
}


#
# Return fusion with the number of samples 
#
countSamples <- function(software,inputDir) { 

  merged.AB <- extractFeatures(software,inputDir)
  
  #############################################################################
  # Summarize the results (gene, #samples, total samples)
  ############################################################################# 
  folders <- list.dirs(path=inputDir,full.names=TRUE,recursive=TRUE)
  numFolders <- length(folders)  
  
  res.AB <- data.frame(geneA=merged.AB$geneA,geneB=merged.AB$geneB,
                       NumSamples=0,TotalSamples=numFolders-1,stringsAsFactors = F)
  res.AB$NumSamples <- rowSums(merged.AB[,-(1:2)])
  
  return(res.AB)
}


# build Decision Tree to separate two groups (e.g., AML and MM)
buildDT <- function(software,inputDirGroup1,inputDirGroup2,outputDir) { 
  feat1 <- extractFeatures(software,inputDirGroup1)
  feat2 <- extractFeatures(software,inputDirGroup2)
  feat <- full_join(feat1,feat2,by=c('geneA','geneB'))
  # replace NA with 0
  feat[is.na(feat)] <- 0
  rownames(feat) <- paste0(feat$geneA,':',feat$geneB)
  feat$geneA <- NULL
  feat$geneB <- NULL
  dat <- data.frame(t(feat),stringsAsFactors = F)
  colnames(dat) <- rownames(feat)
  dat[] <- lapply(dat,as.character)
  dat[] <- ifelse(dat=='1','Y','N')
  label <- gsub('(.+)-.+','\\1',rownames(dat))
  dat$label <- as.factor(label)
  
  # run decision tree and save the result
  fit <- rpart(label ~ ., data=dat, method="class")
  
  outputFilename <- paste0(outputDir,'/DecisionTree.png')
  png(filename=outputFilename)
  prp(fit,varlen=0)
  dev.off()
  cat('Decision tree output saved to',outputFilename,'\n')
}


FisherTest <- function(software,inputDirGroup1,inputDirGroup2,outputDir) { 
  folders1 <- list.dirs(path=inputDirGroup1,full.names=TRUE,recursive=TRUE)
  folders2 <- list.dirs(path=inputDirGroup2,full.names=TRUE,recursive=TRUE)

  if (length(folders1) < 1) {
    stop("The input folder (first group) must contain at least one sample!!")
  } else if (length(folders2) < 1) {
    stop("The input folder (second group) must contain at least one sample!!")
  }
  

  group1 <- countSamples(software,inputDirGroup1)
  Num.group1 <- paste0("Num.",basename(inputDirGroup1))
  Total.group1 <- paste0("Total.",basename(inputDirGroup1))
  colnames(group1) <- c("geneA","geneB",Num.group1,Total.group1)

  group2 <- countSamples(software,inputDirGroup2)
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
  
  
  
  tmp.res.AB <- res %>% select(geneA,geneB)
  # lookup geneA and geneB in other databases
  DB <- c('Human_TSG', 'PanCancer_TSG', 'Oncogene', 'COSMIC', 'Kinase',
          'ChimerKB', 'ChimerPub', 'ChimeSeq'
  )
  for (db in DB) {
    tmp.res.AB <- lookup.DB(tmp.res.AB,db)
  }
  res <- inner_join(res,tmp.res.AB, by=c("geneA","geneB"))  
  
  
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

if (length(args)!=4) {
  stop(paste("\n\n--------------------------------------------------------------------------------------\n\n",
             "Four arguments must be provided.\n",
             "Usage: Rscript co-fuse2.R _SOFTWARE_  _GROUP_1_INPUT_FOLDER_  _GROUP_2_INPUT_FOLDER_  _OUTPUT_FOLDER_ \n\n",
             "where\n",
             "      _SOFTWARE_ is one of 'defuse', 'fusioncatcher', 'tophat', 'soapfuse' or 'generic'\n",
             "      _GROUP_1_INPUT_FOLDER_  is the folder containing the output of fusion software (GROUP 1)\n",
             "      _GROUP_2_INPUT_FOLDER_  is the folder containing the output of fusion software (GROUP 2)\n",
             "      _OUTPUT_FOLDER_ is the folder containing the output of Co-Fuse\n",
             "\n\n--------------------------------------------------------------------------------------\n\n")
       ,call.=FALSE)  
} else {
  cat("Evaluating Fisher's Exact test\n")
  FisherTest(args[1],args[2],args[3],args[4])
  cat("Building Decision tree\n")
  buildDT(args[1],args[2],args[3],args[4])
  cat("------------------------------\n")
}



