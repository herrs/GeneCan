library(tidyverse)
library(zoo)
library(data.table)


getAutoValues <- function(timepoint, GOIdSNP, autodSNP){
  ##create blank dataframe
  timeDF <- data.frame(CHROM = character(), POS = integer() )
  
  ##grabs timepoints from dataframes
  GOIdSNP <- GOIdSNP[grep(timepoint, names(GOIdSNP))]
  autodSNP <- autodSNP[grep(timepoint, names(autodSNP))]
  
  ##create dataframe filled with the position of each SNP from the GOI qtlseqr outputs 
  ## for a single timepoint
  len <- length(GOIdSNP)
  for(i in c(1:len)){
    df <- data.frame(CHROM = GOIdSNP[[i]]$CHROM, POS = GOIdSNP[[i]]$POS)
    timeDF <- rbind(timedf,DF)
  }
  ##remove duplications
  timeDF <- distinct(timeDF)
  
  ##Retrieve location of SNPs and dSNPS for autofluorescence
  aut0 <- data.frame(CHROM =autodSNP[[1]]$CHROM, POS =autodSNP[[1]]$POS, tricubeDeltaSNP =  autodSNP[[1]]$tricubeDeltaSNP)
  ###DF filled with rows not contained in the autofluorescence 
  diff <- setdiff(timeDF, aut0[,1:2])
  
  ##Create nested dataframes for both the GOI and Autofluorescence
  aut0 <- aut0 %>% group_by(CHROM) %>% nest
  diff <- diff %>% group_by(CHROM) %>% nest
  chrom <- (unlist(aut0[,1]))
  GOI <- diff %>% filter(CHROM == chrom[1])
  
  
  ##Call interpolate to estimate dSNP values for snps that are found in GOI but not in autofluorescence for each timepoint
  ##iterates through each chromosome
  z <- lapply(chrom, interpolate, aut0, diff)
  full <- rbindlist(z)
  full$combine <- paste(full$CHROM, full$POS, sep ="_")

  return(full)
  
}


interpolate <- function(chrom, aut0, diff){
  
  #filter on chromsome
  GOI <- diff %>% filter(CHROM == chrom)
  auto <- aut0 %>% filter(CHROM == chrom)
  auto <- auto$data[[1]]
  GOI <- GOI$data[[1]]
  GOI <- GOI[order(GOI$POS),]
  auto <- auto[order(auto$POS),]
  
  ##approximate dSNP value for missing GOI Positions
  if (length(unlist(GOI$POS)) >1){
    imputateAF <- approx(unlist(auto$POS), unlist(auto$tricubeDeltaSNP), xout = unlist(GOI$POS),  rule = 2)
    imputateAF <-  tibble(POS = imputateAF$x, tricubeDeltaSNP = imputateAF$y)
    t <- rbind(auto, imputateAF)
    imputatedValues <- arrange(t, POS)
    ###Replace beginning and end values with nearest value
    imputatedValues$tricubeDeltaSNP <- na.locf(imputatedValues$tricubeDeltaSNP, fromLast =  T, na.rm = F)
    imputatedValues$tricubeDeltaSNP <- na.locf(imputatedValues$tricubeDeltaSNP,  na.rm = F)
    imputatedValues <- cbind(imputatedValues, CHROM = replicate(length(imputatedValues$tricubeDeltaSNP),chrom))
  
  return(imputatedValues)
  }else{
    return(auto)
  }
  
}



