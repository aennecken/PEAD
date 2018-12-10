# written by Philipp Jurmeister (philipp.jurmeister@charite.de) and Anne Schoeler (anne.schoeler@charite.de)
# August 2018

# load required packages
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(gridExtra)
library(randomForest)
library(ggpubr)
library(rmarkdown)
library(Harman)
library(gridExtra)
library(reshape)
library(pander)
library(dplyr)

# Please define your working directory!
# This directory has to include:
# 1. sample_annotation.csv
# 2. Report.Rmd
# 3. all raw IDAT files required for your analysis
# 4. Tree_and_selected_CpGs.Rdata
# The report file will also be saved to this directory

setwd("~/Downloads/PEAD-master/")

# no further adjustments needed after this step, just run the rest of the code
output_dir <- getwd()

# load Random Forest Formula and Probe Names of Probes used as Input for Classification
load("PEAD_metadata.Rdata")

# read in data (specify the directory, where the idat files are located)
targets <- read.metharray.sheet(output_dir, pattern="sample_annotation.csv")

# generate output
# depending on the number of samples the data is either normalized with preprocessSWAN (for <3 samples) or preprocessFunnorm (for 3 or more samples)

if(!exists("cols")&!exists("tree")&!exists("means"))  {
    message("Please load PEAD_metadata.Rdata file")
  
} else {
  if(nrow(targets)<3) {
    # check if Random Forest Formula and Probe Names are loaded 
    # read in raw intensity signals
    rgSet <- read.metharray.exp(targets=targets, force=TRUE)
    
    # normalise data
    mSetSq<- preprocessSWAN(rgSet)
    
    # calculate M-values
    bVals <- getBeta(mSetSq)
    bVals <- shiftBetas(bVals, shiftBy=1e-4)
    mVals <- as.data.frame(logit2(bVals))
    
    # annotate mean values to dataset
    mVals$cpgs <- rownames(mVals)
    mVals.select <- right_join(mVals, means, by="cpgs")
    cpgs <- mVals.select$cpgs
    mVals.select$cpgs <- NULL
    mVals.select <- data.matrix(mVals.select)
    
    # load function to replace NAs with mean values
    remNA <- function(x){
      if(sum(is.na(x)) > 0) x[is.na(x)] <- x[length(x)]
      x
    }
    
    # apply function to dataset
    mVals.select <- as.data.frame(t(apply(mVals.select[,],1,remNA)))
    rownames(mVals.select) <- cpgs
    mVals.select$means <- NULL
    mVals.select$cpgs <- NULL
    
    # apply random forest to sample
    Prediction <- predict(tree, t(mVals.select), proximity = TRUE)
    Prediction.prob <- as.data.frame(predict(tree, t(mVals.select), type="prob"))
    
    #generate report
    render("Report.Rmd",output_dir=output_dir)
    message(paste("Your analyis was successfull! The preprocessSWAN function was used to normalise your data. You can find the HTML file containing the results here: "),getwd(),"/Report.html")
    
    
  } else{
    # read in raw intensity signals
    rgSet <- read.metharray.exp(targets=targets, force=TRUE)
    
    # normalise data
    mSetSq<- preprocessFunnorm(rgSet, verbose=FALSE)
    
    # calculate M-values
    bVals <- getBeta(mSetSq)
    bVals <- shiftBetas(bVals, shiftBy=1e-4)
    mVals <- as.data.frame(logit2(bVals))
    
    # annotate mean values to dataset
    mVals$cpgs <- rownames(mVals)
    mVals.select <- right_join(mVals, means, by="cpgs")
    cpgs <- mVals.select$cpgs
    mVals.select$cpgs <- NULL
    mVals.select <- data.matrix(mVals.select)

    # load function to replace NAs with mean values
    remNA <- function(x){
      if(sum(is.na(x)) > 0) x[is.na(x)] <- x[length(x)]
      x
    }
    
    # apply function to dataset
    mVals.select <- as.data.frame(t(apply(mVals.select[,],1,remNA)))
    rownames(mVals.select) <- cpgs
    mVals.select$means <- NULL
    mVals.select$cpgs <- NULL

    # apply random forest to sample
    Prediction <- predict(tree, t(mVals.select), proximity = TRUE)
    Prediction.prob <- as.data.frame(predict(tree, t(mVals.select), type="prob"))
    
    # generate report
    render("Report.Rmd",output_dir=output_dir)
    message(paste("Your analyis was successfull! The preprocessFunnorm function was used to normalise your data. You can find the HTMIL file containing the results here: "),getwd(),"/Report.html")
   }
 }