library(tidyverse)
library(jsonlite)
library(openxlsx)
library(stringr)
















####MAIN FUNCTION###############################################################
###################CHANGE CODE HERE

###Function outputs excel sheets containing possible candidate genes for each hotspot
###ADD gene name after 2nd "_" the gene file name.
###For example for ARO4 the file name is "GOIS/R_ready_AR04_t0_vcf.table
autoQTLsnps <- function(name, snp, genes, region){
  geneName <- unlist(strsplit(name, "_"))[3]
  print(geneName)
  snps <- snp[name]
  region <- region[name]
  qtlDF <- QTLSNPs(snps,genes, region, geneName)
  new <- createWorkbook()
  nam <- names(qtlDF)
  
  Map(function(data, nam){
    addWorksheet(new, nam)
    
    writeData(new, nam, data)
    
  }, qtlDF, nam)
  
  saveWorkbook(new, file = paste(substring(name, 14, 25), ".xlsx", sep =""), overwrite = TRUE)
  return(qtlDF)
  
}


##Function inputs location of SNP list, location of Genes and location of sig. QTLs
##Function outputs genes within the significant region with the average delta SNP value
##INPUT for SNP must be a folder container the different experiments
QTLSNPs<- function(SNP, Gene, QTL, GOI){
  
  QTLs <- list()
  snpFiles <- SNP[[1]]
  
  sig_Regions <- QTL[[1]]
  sig_Regions$CHROM <- as.character(sig_Regions$CHROM)
  
  
  len <- lengths(sig_Regions)[1]
  Genes <- read.delim(Gene,header=FALSE)
  Genes <- getGenes(Genes, sig_Regions)
  chromLen <- getChromLen(Genes)
  Genes <- subset(Genes,Region=="gene")
  geneList <- geneListF(Genes,sig_Regions)
  
  
  dSNPList <- getSNP_INA1(snpFiles, sig_Regions)
  for(i in 1:len){
    ##retrieve chrNum of QTL region
    chrNum <- sig_Regions[i,1]
    ##get the start and end of qtl to be used in f2
    QTLstart <- sig_Regions[i,3] 
    QTLend <- sig_Regions[i,4]
    ##get the SNPs from the QTL region
    chromSNP <- as.data.frame(getChr(dSNPList, chrNum), col.names = NULL)
    
    ##Get variables to use in the f1 function
    ##where chromosome starts, ends and the last row of gene table
    chromlength = chromLen %>% filter(CHROM == chrNum)
    chromstart <- chromlength[2]
    chromend <- chromlength[3]
    lastrow <- lengths(geneList[[chrNum]])[1]
    chrom <- as.data.frame(geneList[names(geneList) == chrNum])
    if(as.numeric((substring(chrNum,4))) < 10){
      colnames(chrom) <- substring(colnames(chrom), 6)
    }else{
      colnames(chrom) <- substring(colnames(chrom), 7)
    }
    
    ##Creat gene regions with f1
    Chrom <- f1(chrom, chromstart, chromend, lastrow)
    ##combine the different experiments, selecting only the needed columns
    
    dSNP <- chromSNP %>% dplyr::select(CHROM,POS, tricubeDeltaSNP)
    SNPnames <- c("CHROM", "POS", "dSNP")
    colnames(dSNP) <- SNPnames
    ##Reduce gene list to only those within the QTL region
    ChromQTL1 <- f2(Chrom, QTLstart, QTLend)
    qLen <- lengths(ChromQTL1)[1]
    if (dim(ChromQTL1)[1] == 0){
      next
    }
    
    
    ChromGeneSNPdf <- f3(dSNP,ChromQTL1)
    ChromGeneGdf <- f4(ChromGeneSNPdf,ChromQTL1)
    ChromGeneGdf <- ChromGeneGdf %>% mutate_at(vars(dSNP), f5)
    QTLs[[i]] <- ChromGeneGdf
  }
  name <- paste(sig_Regions$CHROM,as.character(sig_Regions$qtl),sep ="_")
  names(QTLs) <- name
  
  
  allGenes <- lapply(QTLs, GeneCan, GOI = GOI)
  unNest <- map(allGenes, unnested)
  unList <- map(unNest, unlisted)
  QTLs <- lapply(QTLs, getChromName)
  name <- names(QTLs)
  print(name)
  len <- length(QTLs)
  df <- list()
  for (i in 1:len){
    df[[i]] <- inner_join(QTLs[[i]], unList[[i]],by = "geneName")
  }
  names(df) <- name
  return(df)
}




#import files
ListFiles <- function(path,pattern = ".csv") {
  require(gtools)
  myFiles <- list.files(path="SNPs", pattern=pattern, full.names=TRUE, recursive=FALSE)
  mixedsort(myFiles)
}


########################code for dSNP portion##########################

##get SNPS
getSNP_INA1 <- function(SNPs, sigRegion){
  SNPsChrom <-  geneListF(SNPs, sigRegion) 
  return(SNPsChrom)
  
}

##For the length of each chromosome
getChromLen <- function(genes) {
  chromLen <- (genes %>% filter(substr(CHROM,1,19) == "##sequence-region C"))[1]
  chromLen <- substr(chromLen$CHROM, 19,33)
  chromLen <- data.frame(do.call(rbind, (strsplit(chromLen, " "))))
  colnames(chromLen ) <- c("CHROM", "start", "end")
  chromLen$start <- as.integer(chromLen$start)
  chromLen$end <- as.integer(chromLen$end)
  return(chromLen)
}

#trim the gene list to only include ORFs for a given chromosome
getGenes <- function(gene, sigRegion){
  gene <- gene %>% dplyr::select(c("V1","V3","V4","V5","V6","V7","V8","V9","V11"))
  names(gene) <- c("CHROM", "Region", "ORFstart", "ORFend", "SoR", "EoR","Range","GeneID","Name")
  gene$Region <- as.character(gene$Region)
  gene$ORFstart <- as.numeric(gene$ORFstart)
  gene$ORFend <- as.numeric(gene$ORFend)
  return(gene)
  
}

#create a subset of the gene list to only include genes for a given chromosome
#create a subset of the SNP list to only include SNPs for a given chromosome--sub in different numbers depending on chromosome
geneListF <- function(Genes, SigR){
  geneList <- list()
  chrom = unique(SigR$CHROM)
  len = length(chrom)

  for (i in 1:len){
    geneList[[i]] <-  tryCatch(
      {subset(Genes, Chromosome == chrom[i])},
      error = function(cond){
        snp <- subset(Genes, CHROM == chrom[i]) 
        return (snp)
      } 
    )
  }
  names(geneList) = chrom
  return(geneList)
}

#STEP 2: CREATE FUNCTION TO MAKE GENE REGIONS
f1 <- function(chromosome, chromstart,chromend, lastrow) {
  new <- chromosome
  new$SoR <- "NA"
  new$EoR <- "NA"
  new$Range <- "NA"
  lenSoR <- length(new$SoR)
  for (r in 1:lenSoR) {
    if (r > 1) {
      new$SoR[r] <- (new$ORFend[r-1])
    }
  }
  new$SoR[1] <- chromstart
  new$SoR <- as.numeric(new$SoR)
  lenEoR <- length(new$EoR)
  for (l in 1:lenEoR) {
    if (l < lastrow) {
      new$EoR[l] <- (new$ORFstart[l+1])
    }
  }
  new$EoR[lastrow] <- chromend
  new$EoR <- as.numeric(new$EoR)
  lenRan <- length(new$Range)
  for (i in 1:lenRan) {
    new$Range[i] <- (new$EoR[i]-new$SoR[i])
  }
  new$Range <- as.numeric(new$Range)
  chromosome$SoR <- new$SoR
  chromosome$EoR <- new$EoR
  chromosome$Range <- new$Range
  return(chromosome)
}

#STEP 3: CREATE FUNCTION TO REDUCE GENE LIST TO JUST QTL AREA
f2 <- function(chromosome,QTLstart, QTLend) {
  QTL <- subset(chromosome,EoR>=QTLstart & SoR<=QTLend)
}
getChr <- function(data, chr){
  chrom <- (data[names(data) == chr])
  return(chrom)
}

#STEP 4: MAKE A FUNCTION TO START SORTING SNP
#create a new dataset that includes gene, gene region, all associated SNPs and G' scores for all of those SNPs

f3 <- function(chromosomeSNPs,chromosomeqtl) {
  GeneSNPdf <- data.frame(geneName = factor(), snp.coord = integer(),"dSNP" = integer())
  snpstart <- 1
  lastrowadded <- integer(length(chromosomeqtl$Name))
  len <- length(chromosomeqtl$Name)
  for (g in 1:len) {
    for (s in snpstart:length(chromosomeSNPs$POS)) {
      if (chromosomeqtl$SoR[g] <= chromosomeSNPs$POS[s] & chromosomeSNPs$POS[s] <= chromosomeqtl$EoR[g]) {
        
        GeneSNPdf <- rbind(GeneSNPdf, data.frame(geneName = chromosomeqtl$Name[g], SoR=chromosomeqtl$SoR[g], EoR=chromosomeqtl$EoR[g], snp.coord=chromosomeSNPs$POS[s],"dSNP" = chromosomeSNPs$dSNP[s]))
      }
      else if (chromosomeqtl$EoR[g] < chromosomeSNPs$POS[s]) {
        lastrowadded[g] <- s
        break
      }
    }
    if (g > 1) {
      if (lastrowadded[g - 1] > 0) {
        snpstart <- lastrowadded[g - 1]
      }
    }
  }
  return(GeneSNPdf)
}
#STEP 5: MAKE A FUNCTION TO CREATE GENE GPRIME DF
##I should make this in fewer lines but its done so...

f4 <- function(GeneSNPdf,chromosomeqtl) {
  m2t2 =aggregate(GeneSNPdf$dSNP,list(GeneSNPdf$gene,GeneSNPdf$SoR,GeneSNPdf$EoR),mean)
  names(m2t2)[names(m2t2)=="Group.1"] <- "geneName"
  names(m2t2)[names(m2t2)=="Group.2"] <- "SoR"
  names(m2t2)[names(m2t2)=="Group.3"] <- "EoR"
  names(m2t2)[names(m2t2)=="x"] <- "dSNP"
  geneG <- m2t2
  newG <- rbind(geneG,data.frame(geneName=chromosomeqtl$Name,SoR=chromosomeqtl$SoR,EoR=chromosomeqtl$EoR,dSNP = "NA"))
  GeneGdf <- newG[!duplicated(newG$gene),]
  GeneGdf <- GeneGdf %>% arrange(SoR)
  GeneGdf$dSNP <- as.numeric(GeneGdf$dSNP)
 
  return(GeneGdf)
}

#for genes with "NA" as the mean G', take the mean of the closest gene G' above and below
f5 <- function(dat) {
  N <- length(dat)
  new <- dat
  na.pos <- which(is.na(dat))
  if (length(na.pos) %in% c(0, N)) {
    return(dat)
  }
  non.na.pos <- which(!is.na(dat))
  for (i in 1:length(na.pos)) {
    left.pos <- max(non.na.pos[non.na.pos < na.pos[i]])
    right.pos <- min(non.na.pos[non.na.pos > na.pos[i]])
    if (left.pos < 1) {left.pos <- right.pos}
    if (right.pos > N) {right.pos <- left.pos}
    new[na.pos[i]] <- mean(c(dat[left.pos], dat[right.pos]))
  }
  return(new)
}

################################Functions for go terms and phenotypes###################
##get chromosome list, enter path to chromosome file
geneTable <- function(chrom){
  return(as.data.frame(read.delim(chrom ,sep = "")))
  
}

##input chromosome list and gene that you want compared to
GeneCan <- function(chromosomeList, GOI){
  
  
  ##Get interactions and regulators from GOI
  
  regGOI <- fromJSON(paste("https://www.yeastgenome.org/backend/locus/" ,GOI, "/regulation_details", sep =""))
  interGOI <-  fromJSON(paste("https://www.yeastgenome.org/backend/locus/" ,GOI, "/interaction_details", sep =""))
  interGOI <- interGOI$locus2$display_name
  regGOI <- regGOI$locus1$display_name
  
  #######Phenotypic and go search terms####################################################
  ###CHANGE THE TERMS HERE######################################################################
  flags <- c("MAPK",
             "mating",
             "pheromone response",
             "cell cycle",
             "protein decay",
             "translation efficiency",
             "translation initiation",
             "ubiquitin mediated",
             "ubiquitination",
             "ubiquitin",
             "polysome",
             "shmoo",
             "translation regulation",
             "ribosome",
             "translation",
             "translate",
             "tRNA",
             "gene expression",
             "posttranscriptional",
             "cell fusion",
             "regulation of cell size",
             "cell size",
             "protein/peptide accumulation",
             "G1")
  
  
  ##Filter the chromsome list to just get the gene names
  geneList <- c()
  chromosomeList <- getChromName(chromosomeList)
  geneList <- unlist(chromosomeList$geneName)
  print(geneList)
  
  ###"Slow" for loop to access every gene to get the information and put them in a tibble
  len = length(geneList)
  if(len > 1){
    for (j in 2:len){
      if ( j > 2){
        geneInfo <- bind_rows(geneInfo, findInfo(geneList[j], interGOI, regGOI, flags))
      }else {
        geneInfo <- bind_rows(findInfo(geneList[(j-1)], interGOI, regGOI,flags), findInfo(geneList[(j)], interGOI, regGOI, flags))
      }
    }
  }else if (len == 1){
    
    geneInfo <- as_tibble(findInfo(geneList[(len)], interGOI, regGOI, flags))
  }else{
    return(NA)}
  
  ###retrieves the total number of matches each gene has
  geneInfo1 <- geneInfo %>% group_by(geneName) %>% nest
  geneInfo2 <- geneInfo1 %>% mutate(count = map(data, getLengths))
  geneInfo2$count <- unlist(geneInfo2$count)
  ###Sorts them based on number
  sortedGenes <- geneInfo2 %>% arrange(desc(count))
  return(sortedGenes)
  
}



###Gathers JSON data from yeastGenome.org
##Finds the genes from a chromsomes interactions with the interactors/regulators with FIG1 
findInfo <- function(geneList, interGOI, regGOI, flags){    
  
  geneName = geneList
  
  pheno <- fromJSON(paste("https://www.yeastgenome.org/backend/locus/" ,geneName, "/phenotype_details", sep =""))
  inter <- fromJSON(paste("https://www.yeastgenome.org/backend/locus/" ,geneName, "/interaction_details", sep =""))
  
  
  reg <- fromJSON(paste("https://www.yeastgenome.org/backend/locus/" ,geneName, "/regulation_details", sep =""))
  
  go <- fromJSON(paste("https://yeastgenome.org/backend/locus/", geneName,"/go_details", sep = ""))
  
  inter <- inter$locus2$display_name
  reg <- reg$locus1$display_name
  pheno <- pheno$phenotype$display_name
  go <- go$go$display_name
  
  goF <- unique(grep(paste(flags,collapse="|"), 
                     go, value=TRUE))
  phenoF <-  unique(grep(paste(flags,collapse="|"), 
                         pheno, value=TRUE))
  
  interIF <- inter[inter %in% interGOI]
  interRF <- inter[inter %in% regGOI]
  interF <- unique(c(interIF, interRF))
  regRF <- reg[reg %in% regGOI]
  regIF <- reg[reg %in% interGOI]
  regF <- c(regRF, regIF)
  regF <- check0(regF)
  interF <- check0(interF)
  goF <- check0(goF)
  phenoF <- check0(phenoF)
  
  geneRow <- tibble(geneName =  geneName, phenotypes = list(phenoF), interactions = list(interF),regulators = list(regF), goterms = list(goF))
  return(geneRow)                 
}  

###This function retrieves the count of matches found within each gene

getLengths <- function(val){
  vec <- c()
  for (i in 1:4){
    vec[i] <- length(unlist(val[i]))
  } 
  return(sum(vec))
}



##This function is to make sure there is a value in every column for each gene
check0 <- function(value){
  if (identical(value, character(0)) | is.null(value)){
    return ("no value")
  }else{
    return(value)
  }
}
##Get chromosome number to name the tables
getChrom <- function(chroms){
  chroms <- substring(chroms, 7,20)
  return(chroms)
  
}

getChromName <- function(chromosomeList){
  n <- chromosomeList$geneName
  chromosomeList$geneName <- substring(n, 6)
  return(chromosomeList)
}

###Formatting functions for exporting to excel sheet

unnested <- function(chromT){
  new <- chromT  %>% unnest(col = "data")
  return(new)
}

unlisted <- function(try){
  try$phenotypes <-  map(try$phenotypes, toString) %>% unlist(try$phenotypes)
  try$interactions <-  map(try$interactions, toString) %>% unlist(try$interactions)
  try$regulators <-  map(try$regulators, toString) %>% unlist(try$regulators)
  try$goterms <-  map(try$goterms, toString) %>% unlist(try$goterms)
  return(try)
}





