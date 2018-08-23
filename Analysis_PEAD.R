# written by Philipp Jurmeister (philipp.jurmeister@charite.de) and Anne Schoeler (anne.schoeler@charite.de)
# August 2018

# load required packages
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gridExtra)
library(randomForest)
library(ggpubr)
library(rmarkdown)
library(Harman)
library(gridExtra)
library(reshape)

# set working directory, this directory has to include:
# 1. sample_annotation.csv
# 2. Report.Rmd
# 3. all raw IDAT files required for your analysis
# 4. Tree_and_selected_CpGs.Rdata
# The report file will also be saved to this directory
setwd("/Volumes/Data/PEAD/")

output_dir <- getwd()

# load Random Forest Formula and Probe Names of Probes used as Input for Classification
load("Tree_and_selected_CpGs.Rdata")

# read in data (specify the directory, where the idat files are located)
targets <- read.metharray.sheet(output_dir, pattern="sample_annotation.csv")

# generate output
# depending on the number of samples the data is either normalized with preprocessSWAN (for <3 samples) or preprocessFunnorm (for 3 or more samples)

if(!exists("cols")&!exists("tree"))  {
    message("Please load Tree_and_selected_CpGs.Rdata file")
  
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
    mVals <- logit2(bVals)
    
    # grab the 10000 Probes from the random forest training
    mVals.select <- as.data.frame(mVals[c(cols),])
    
    # apply random forest to sample
    Prediction <- predict(tree, t(mVals.select), proximity = TRUE)
    Prediction.prob <- predict(tree, t(mVals.select), type="prob")
    
    #generate report
    render("Report.Rmd",output_dir=output_dir)
    message(paste("Your analyis is finished. The preprocessSWAN function was used to normalise your data. You can find the PDF containing the results here: "),getwd(),"/Report.html")
    
    
  } else{
    # read in raw intensity signals
    rgSet <- read.metharray.exp(targets=targets, force=TRUE)
    
    # normalise data
    mSetSq<- preprocessFunnorm(rgSet, verbose=FALSE)
    
    # calculate M-values
    bVals <- getBeta(mSetSq)
    bVals <- shiftBetas(bVals, shiftBy=1e-4)
    mVals <- logit2(bVals)
    
    # grab the 10000 Probes from the random forest training
    mVals.select <- as.data.frame(mVals[c(cols),])
    
    # apply random forest to sample
    Prediction <- predict(tree, t(mVals.select), proximity = TRUE)
    Prediction.prob <- as.data.frame(predict(tree, t(mVals.select), type="prob"))
    
    # generate report
    render("Report.Rmd",output_dir=output_dir)
    message(paste("Your analyis is finished. The preprocessFunnorm function was used to normalise your data. You can find the PDF containing the results here: "),getwd(),"/Report.html")
   }
 }







