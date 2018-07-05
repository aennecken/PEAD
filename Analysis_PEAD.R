# written by Philipp Jurmeister (philipp.jurmeister@charite.de) and Anne Schoeler (anne.schoeler@charite.de)
# July 2018

# load required packages
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gridExtra)
library(randomForest)
library(Harman)

# set working directory, the output PDF is saved in the same folder
setwd("/Users/admin/Documents/Kollaborationen/Philipp Jurmeister/PEAD")

# load Random Forest Formula and Probe Names of Probes used as Input for Classification
load("Tree_and_selected_CpGs.Rdata")

# read in data (specify the directory, where the idat files are located)
targets <- read.metharray.sheet("/Users/admin/Documents/Kollaborationen/Philipp Jurmeister/PEAD", pattern="sample_annotation.csv")

# generate output
# depending on the number of samples the data is either normalized with preprocessSWAN (for <3 samples) or preprocessFunnorm (for 3 or more samples)

if(!exists("cols")&!exists("tree"))  {
    message("Please load Tree_and_selected_CpGs.Rdata file")
  
} else {
  if(nrow(targets)<3) {
    # check if Random Forest Formula and Probe Names are loaded 
    # read in raw intensity signals
    rgSet <- read.metharray.exp(targets=targets)
    
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
    
    idat_name<-gsub(".*/","",targets$Basename)
    output <- data.frame(Sample_ID=targets$sampleID,idat=idat_name,Random_forest_prediction=Prediction$predicted, Probability_CRAD=Prediction.prob[,1],Probability_LUAD=Prediction.prob[,2])
    rownames(output)<-NULL
    pdf("Random_forest_Output_Swan.pdf", height=8.5, width=11)
    grid.table(output, rows = NULL)
    dev.off()
    message(paste("Your analyis is finished. The preprocessSWAN function was used to normalise your data. You can find the PDF containing the results here:"),getwd(),"/Random_forest_Output_Swan.pdf")
    
    
  } else{
    # read in raw intensity signals
    rgSet <- read.metharray.exp(targets=targets)
    
    # normalise data
    mSetSq<- preprocessFunnorm(rgSet, verbose = FALSE)
    
    # calculate M-values
    bVals <- getBeta(mSetSq)
    bVals <- shiftBetas(bVals, shiftBy=1e-4)
    mVals <- logit2(bVals)
    
    # grab the 10000 Probes from the random forest training
    mVals.select <- as.data.frame(mVals[c(cols),])
    
    # apply random forest to sample
    Prediction <- predict(tree, t(mVals.select), proximity = TRUE)
    Prediction.prob <- predict(tree, t(mVals.select), type="prob")
    
    idat_name<-gsub(".*/","",targets$Basename)
    output <- data.frame(Sample_ID=targets$sampleID,idat=idat_name,Random_forest_prediction=Prediction$predicted, Probability_CRAD=Prediction.prob[,1],Probability_LUAD=Prediction.prob[,2])
    rownames(output)<-NULL
    pdf("Random_forest_Output_Funnorm.pdf", height=8.5, width=11)
    grid.table(output, rows = NULL)
    dev.off()
    message(paste("Your analyis is finished. The preprocessFunnorm function was used to normalise your data. You can find the PDF containing the results here:"),getwd(),"/Random_forest_Output_Funnorm.pdf")
   }
 }







