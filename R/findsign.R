library(tidyverse)

normSigR[[17]] <- normSigRfig[[4]]
normSigR[[18]] <- normSigRfig[[8]]
names(normSigR)[17] <- names(normSigRfig)[4]
names(normSigR)[18] <- names(normSigRfig)[8]


iterate <- names(normSigR)
negative <- list()
positive <- list()


splitNeg <- function(name, data){
  try <- data[[name]]
  neg <- try %>% filter(avgDeltaSNP <  0)
  return(neg)
}

splitPos <- function(name, data, positive){
  
  try <- data[[name]]
  pos <- try %>% filter(avgDeltaSNP > 0)
  return(pos)
}

negative <- lapply(iterate, splitNeg, normSigR)
positive <- lapply(iterate, splitPos, normSigR)
names(negative) <- iterate
names(positive) <- iterate
