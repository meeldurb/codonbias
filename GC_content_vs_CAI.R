#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the mean CAI vs the 
##GC content of all genomes
###################################################################

# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")


# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

genomecount = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  domain.files <- paste("Domain_data_ENA/", genomeID, ".csv", sep = "")
  if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)
    seq.data <- read.csv(file = domain.files, sep = ",", header = TRUE, as.is = TRUE)
    
    # calculate mean of all the CAI values in the genome
    mean.cai <- mean(cai.data[,2])
    
    # calculate the GC content of the genome
    all.seq <- paste(as.matrix(seq.data)[,4], sep="", collapse="")
    seq.split <- strsplit(all.seq, "")[[1]]
    GCcont <- GC(seq.split)*100
    xlim = c(10, 90)
    ylim = c(0.1, 1)
    # then plot the GC content against the mean CAI
    if (genomecount == 0){
      plot(GCcont, mean.cai, type = "p", xlim = xlim, ylim = ylim,
           pch = 18, col = "blue",
           main = "Average CAI vs. GC content",
           xlab = "GC content (%)", 
           ylab = "Mean CAI")
      grid(NULL, NULL, lty = 6, col = "cornsilk2")
      genomecount = genomecount + 1
    } 
    else {
      points(GCcont, mean.cai, pch = 18, col = "blue")
    } 
  }
}



######______________________________For 1 genome______________________________######

genomeID <- "GCA_000003645"
cai.files <- paste("CAI_complete_ENA/", genomeID, "_CAI_complete.csv", sep = "")
cai.data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)
domain.files <- paste("Domain_data_ENA/", genomeID, ".csv", sep = "")
seq.data <- read.csv(file = domain.files, sep = ",", header = TRUE, as.is = TRUE)

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


mean.cai <- mean(cai.data[,2])


all.seq <- paste(as.matrix(seq.data)[,4], sep="", collapse="")
GC(all.seq)
seq.split <- strsplit(all.seq, "")[[1]]
GCcont <- GC(seq.split)*100
#nfreq <- table(seq.split)
#GCcont <- ((nfreq[2] + nfreq[3])/sum(nfreq))*100

