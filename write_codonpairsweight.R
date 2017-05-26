#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for drawing data of all the weight values
##of all genomes to a dataframe
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000008185"
genomeID <- "GCA_000008205"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genome.and.organisms <- genome.and.organisms[,1:2]



####______________ create dataframe of all genomes with the codonpairs and their weights of the pairs ______________####
genomecount = 0
n = 0
c = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")  
  w.files <- paste("codonpairs_itweight/", genomeID, "_witcodpairs.csv", sep = "") 
  if (file.exists(w.files)){
    # Searching for Mycoplasma and Spiroplasma, which use other genetic codon table
    match.words <- c("Mycoplasma", "Spiroplasma")
    # i contains the indices where Myco/Spiro is found
    i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
    genomesMycoSpiro <- genome.and.organisms[i,1]
    # when the genome number of myc/spir is not found it will combine the genomeID in the dataframe
    # else the dataframe is too short (less stop codons are omitted)
    # so it will not complain anymore
    if(!(genomeID %in% genomesMycoSpiro)) {
      n <- n + 1
      w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
      # order on codon to get better readability
      ordered.w <- w.data[with(w.data, order(w.data[,1])), ]
      # removing all codon pairs which start with stop codon
      # for weight tables
      stop.pair.i <- grep("\\*_", w.data[,1])
      length(w.data[,1]) - length(stop.pair.i)
      w.pair <- w.data[-stop.pair.i, ]
      length(w.pair[,1])
      # only leaving numbers
      w <- w.pair[,3]
      length(w)
      # getting dataframe of codons and weight values
      if (genomecount == 0){
        codgendf <- data.frame(row.names=w.pair[,2])
        codgendf <- cbind(codgendf, w)
        colnames(codgendf)[n] <- genomeID
        genomecount = genomecount + 1
      } else {
        codgendf <- cbind(codgendf, w)
        colnames(codgendf)[n] <- genomeID
      } 
    } else { 
      print(paste("genomefile", genomeID, "is myco/spiro genome"))
      # c <- c + 1
    }
  } else {
    print(paste("genomefile", genomeID, "does not exist"))
    #  c <- c + 1
  }
}
save(codgendf, file = "codonpairsGenomeDataSet.RData")