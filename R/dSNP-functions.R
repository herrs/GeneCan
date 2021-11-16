

#load the package
library("QTLseqr")
library("tidyverse")
library("gtools")

source("autoNormal.R")



limitIt <- function(df){
  df$tricubeDeltaSNP <- -1* df$tricubeDeltaSNP
  return(df)
}

ListFiles <- function(path,pattern = ".table") {
  require(gtools)
  myFiles <- list.files(path=path, pattern=pattern, full.names=TRUE, recursive=FALSE)
  mixedsort(myFiles)
}

getTime <- function(time, GOI, AF){
  timeGOI <- GOI[grepl(time,GOI)]
  timeAF <- AF[time]
  x <- lapply(timeGOI, AutoCor,  timeAF)
  names(x) <- timeGOI
  return(x)
}


###GET autofluorescence DSNP scores
AutoDSNP <- function(auto){
  #Set sample and file names
  HighBulk <- "High"
  LowBulk <- "Low"
  
  #Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
  Chroms <- paste0(rep("Chr", 17), 1:16)
  
  #Import SNP data from file
  af <-
    importFromGATK(
      file = auto,
      highBulk = HighBulk,
      lowBulk = LowBulk,
      chromList = Chroms
    )
  
 
  
  ##view data to assess filtering options
  minD <- quantile(af$DP.HIGH, probs = seq(0,1,0.05), na.rm = T)[2] + quantile(af$DP.LOW, probs = seq(0,1,0.05), na.rm =T)[2]
  maxD <- quantile(af$DP.HIGH, probs = seq(0,1,0.05), na.rm =T)[20] + quantile(af$DP.LOW, probs = seq(0,1,0.05), na.rm =T)[20]
  
  
  #Filter SNPs based on some criteria
  af_data_filt <-
    filterSNPs(
      SNPset = af,
      refAlleleFreq = 0.2,
      minTotalDepth = minD,
      maxTotalDepth = maxD,
      minSampleDepth = 100,
      minGQ = 99
    )


  
  af_QTL_filt <- runQTLseqAnalysis(
    af_data_filt,
    windowSize = 3e4,
    popStruc = "F2",
    bulkSize = 5000,
    intervals = c(90,95, 99))
  
  af_QTL_filt$combine <- paste(af_QTL_filt$CHROM, af_QTL_filt$POS, sep ="_")
  return(af_QTL_filt)
}
###GET GOI DSNP scores
GOIdSNP <- function( goi){
  
  #Set sample and file names
  HighBulk <- "High"
  LowBulk <- "Low"
  
  #Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
  Chroms <- paste0(rep("Chr", 17), 1:16)
  
  #Import SNP data from file
  
  GOI <- importFromGATK(
    file = goi,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )
  
  
  ##view data to assess filtering options
  
  
  #Filter SNPs based on some criteria
  minD <- quantile(GOI$DP.HIGH, probs = seq(0,1,0.05), na.rm =TRUE)[2] + quantile(GOI$DP.LOW, probs = seq(0,1,0.05), na.rm =TRUE)[2]
  maxD <- quantile(GOI$DP.HIGH, probs = seq(0,1,0.05), na.rm =TRUE)[20] + quantile(GOI$DP.LOW, probs = seq(0,1,0.05), na.rm =TRUE)[20]
  
  goi_data_filt <-
    filterSNPs(
      SNPset = GOI,
      refAlleleFreq = 0.2,
      minTotalDepth = minD,
      maxTotalDepth = maxD,
      minSampleDepth = 100,
      minGQ = 99
    )
  
  
  
  
  
  
  GOI_QTL_filt <- runQTLseqAnalysis(goi_data_filt,
                                    windowSize = 3e4,
                                    popStruc = "F2",
                                    bulkSize = 5000,
                                    intervals = c(90, 95, 99))
  
  return(GOI_QTL_filt)
}



#########Auto fluorescence occurs here
AutoCor <- function( goi, af){
  af_QTL_filt <- af[[1]]
  #Set sample and file names
  HighBulk <- "High"
  LowBulk <- "Low"

  #Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
  Chroms <- paste0(rep("Chr", 17), 1:16)
  
  #Import SNP data from file
  
  GOI <- importFromGATK(
    file = goi,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )
  
  
  ##view data to assess filtering options
  
  
  #Filter SNPs based on some criteria
  minD <- quantile(GOI$DP.HIGH, probs = seq(0,1,0.05), na.rm = T)[2] + quantile(GOI$DP.LOW, probs = seq(0,1,0.05), na.rm =T)[2]
  maxD <- quantile(GOI$DP.HIGH, probs = seq(0,1,0.05), na.rm =T)[20] + quantile(GOI$DP.LOW, probs = seq(0,1,0.05), na.rm =T)[20]
  
  goi_data_filt <-
    filterSNPs(
      SNPset = GOI,
      refAlleleFreq = 0.2,
      minTotalDepth = minD,
      maxTotalDepth = maxD,
      minSampleDepth = 100,
      minGQ = 99
    )
  
  
  
 
  
  GOI_QTL_filt <- runQTLseqAnalysis(goi_data_filt,
                                    windowSize = 3e4,
                                    popStruc = "F2",
                                    bulkSize = 5000,
                                    intervals = c(90, 95, 99))
  
  
  GOI_QTL_filt$combine <- paste(GOI_QTL_filt$CHROM, GOI_QTL_filt$POS, sep ="_")
  x <- within(merge(GOI_QTL_filt,af_QTL_filt,by="combine"), {
    newcube <- tricubeDeltaSNP.x - tricubeDeltaSNP.y
  })
  
  GOI_QTL_filt_3 <- left_join(GOI_QTL_filt, x)
  GOI_QTL_filt_3$tricubeDeltaSNP <- coalesce(GOI_QTL_filt_3$newcube,GOI_QTL_filt_3$tricubeDeltaSNP)
  return(GOI_QTL_filt_3)
}





##see how filtering affected above distributions

  
getQTLTables <-
   function(SNPset,
            method = "QTLseq",
            alpha = 0.05,
            interval = 95){
   
     #QTL <- getSigRegions(SNPset = SNPset, method = method, interval = interval, alpha = alpha)
     SNPset <- SNPset %>%
       dplyr::group_by(CHROM)
     
     if (method == "QTLseq") {
       qtltable <-
         SNPset %>% 
         dplyr::mutate(passThresh = abs(tricubeDeltaSNP) > abs(CI_95)) %>%
         dplyr::group_by(CHROM, run = {
           run = rle(passThresh)
           rep(seq_along(run$lengths), run$lengths)
         }) %>%
         dplyr::filter(passThresh == TRUE) %>% 
         dplyr::ungroup() %>%
         dplyr::group_by(CHROM) %>% 
         dplyr::group_by(CHROM, qtl = {
           qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
         }) %>%
         #dont need run variable anymore
         dplyr::select(-run) %>%
         dplyr::summarize(
           start = min(POS),
           end = max(POS),
           length = end - start,
           nSNPs = length(POS),
           avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
           peakDeltaSNP = ifelse(
             mean(tricubeDeltaSNP) >= 0,
             max(tricubeDeltaSNP),
             min(tricubeDeltaSNP)
           ),
           posPeakDeltaSNP = POS[which.max(abs(tricubeDeltaSNP))],
           avgDeltaSNP = mean(tricubeDeltaSNP)
         )
     } 
     
     qtltable <- as.data.frame(qtltable)
     
  
     return(qtltable)
     
     
   }
 
 
getImputation <- function(allGOI, allAf, times){
  dSNPs <- list()
  dSNPS <- lapply(allGOI,GOIdSNP)
  names(dSNPS) <- allGOI
  autodSNP <- lapply(allAf, AutoDSNP)
  names(autodSNP) <-  times
  imputatedAuto <- lapply(times, getAutoValues, dSNPS, autodSNP)
  names(imputatedAuto) <- c( "t0", "t2","t5", "t8")
  return(imputatedAuto)
  
}




