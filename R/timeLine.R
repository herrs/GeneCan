library(tidyverse)
library(locfit)
library(data.table)
library(rlist)



#plot slope and CI

t0<- read.csv("Filtered_Gprime_SNPs_af_t0.csv")
t2<- read.csv("Filtered_Gprime_SNPs_af_t2.csv")
t5<- read.csv("Filtered_Gprime_SNPs_af_t5.csv")
t8<- read.csv("Filtered_Gprime_SNPs_af_t8.csv")
 
auto <- list(t0,t2,t5,t8)

auto <-auto %>% map( ~select(.x, c(CHROM,POS, tricubeDeltaSNP)))
 
autoFull <- auto %>% purrr::reduce( dplyr::full_join, by = c("CHROM","POS"))
names(autoFull)[3:6] <- c("0","2","5","8" )
autoFull$CHROM <- as.factor(autoFull$CHROM)
autoFull$POS <- as.factor(autoFull$POS)
autoLong <- autoFull %>% pivot_longer(c("0","2","5","8"), names_to = "time", values_to = "tricubedSNP")
autoLong$time <- as.integer(autoLong$time)
autoLong <-  autoLong %>% group_by(across(c(CHROM,POS))) %>% nest()

mod_fun <- function(df){
  x <- tidy(lm(df$tricubedSNP ~ df$time))
  x$p.value[2]
}

autoModel <- autoLong %>% mutate(model = map(data, mod_fun))
autoModel$model <- unlist(autoModel$model)
summary(autoModel$model)
ggplot(autoModel, aes(y=model)) + geom_boxplot()

################################################################

nam <- names(GOIdSNPS)
namedDF <- lapply(nam, giveName, normSigR)
names(namedDF) <- names(GOIdSNPS)

giveName <- function(Name, df){
  df1 <- df[[Name]]
  x <- cbind(df1, name = replicate(lengths(df1)[1], substr(Name,14,20))) 
  x <- cbind(x, chrom = x$CHROM)
  return(x)
}

         

len <- length(namedDF)
fdf <- data.frame()
for (i in 1:len){
  fdf <- rbind(fdf, namedDF[[i]])
}


sepFunc <- function(key,df){
  new <- df[grepl(key, df$name),]
  new2 <- new %>% group_by(CHROM) %>% nest()
  new2 <- new2[order(fdf2$CHROM),]
  return(new2)
}


key <- c("A", "D", "G","I","X", "t0","t2","t5","t8")

sepList <- lapply(key,sepFunc,fdf)
names(sepList) <- key

getData <- function(p){
  yup <- lapply(p$data, getIntervals)
  return(yup)
}

getIntervals <- function(p){
  if(!is.null(p)){
    count = 1
    len <- lengths(p)[1]
    matches <- list()
    chrom = 1
    chromV <- c()
    for (i in 1:len){
      regions <- data.table::between(p$start[i], p$start, p$end)
      if( sum(regions) > 1){
        regions <- p[regions,]
        if(count >2) {
          if ((sum(matches %in% list(regions)) == 0)){
            print(matches %in% list(regions))
            chromV <- c(chromV, paste(regions$name, collapse =" "))
            matches[[count]] <- regions
            count <- count +1
          }else{
            print(1)
            next
          }
        }else{
          chromV <- c(chromV, paste(regions$name, collapse =" "))
          matches[[count]] <- regions
          count <- count +1
        }
      }
      
    }
    names(matches) <- chromV
    return(matches)
    return()
  }else{
    return()
    }
}

new<- lapply(sepList, getData)
try <- new[[1]]





counter = 0
t4 <- data.frame(CHROM = character(), qtl = integer(), Ar_t0_ = integer(), Ar_t2_ = integer(), Ar_t5_ = integer(), Ar_t8_ = integer(),D_t2_v = integer(), D_t5_v = integer(), D_t8_v = integer(), G_t0_v = integer(), G_t2_v = integer(), G_t5_v = integer(), G_t8_v = integer(), I_5_t0 = integer(), I_m_t0 = integer(), I_m_t2 =integer(), X_t5_v = integer(), X_t8_v = integer() )
t3
t2

firstNest <- function(data,df ){
  x <- lapply(data, nestApply,df )
  
}

nestApply <- function(df,overlap){
 counter <<- counter +1
  overlap <- lapply( df ,convertRegion, overlap, counter)
  return(overlap)
  
}
convertRegion <- function(df2, overlap, Chr){
  name <- df2$name
  try <- df2$qtl
  names(try) <- name
  try1 <- t(as.data.frame(try))
  try1 <- cbind.data.frame(try1, CHROM = paste("Chr", Chr, sep = ""))
  overlap <- full_join(overlap,try1) 
  return(overlap)
}
please <- sapply(new, nestApply, overlap)

please2 <- lapply(please,list.rbind)

please3 <- list.rbind(please2)


library(plyr)
library(rlist)
library(tidyverse)





findAll1 <- function(name, df, cont){
  
  df2 <- df[[name]]
  len <- as.integer(unlist(cont[name]))
  print(len)
  fi <- lapply(df2,findAll2,len)

  return(fi)
  
}

  


findAll2 <- function(df, leng){
 x <- df[which(sapply(df, nrow) == leng)]
 return(x)
}

lenGene <-  list(4,3,4,3,2,4,4,4,4)
names(lenGene) <- key

lenGene[key[1]]


why <- lapply(key, findAll1, new, lenGene)
why2 <- lapply(why, bind_rows)

why3 <- lapply(why2, drop_na)
