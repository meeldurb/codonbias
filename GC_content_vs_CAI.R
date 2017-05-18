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
install.packages("locfit", repos="http://cran.rstudio.com/")


# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")
library("locfit")
# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


genomeIDcol <- NULL
meancaicol <- NULL
GCcontcol <- NULL
#pdf("GCvsCAI_plot_itcount.pdf")
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  CDS.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
  if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, sep = ",", header = TRUE, as.is = TRUE)
    seq.data <- read.csv(file = CDS.files, sep = ",", header = TRUE, as.is = TRUE)
    
    # calculate mean of all the CAI values in the genome
    mean.cai <- mean(cai.data[,2])
    
    # calculate the GC content of the genome
    all.seq <- paste(as.matrix(seq.data)[,2], sep="", collapse="")
    seq.split <- strsplit(all.seq, "")[[1]]
    GCcont <- GC(seq.split)*100
    
    # fill columns with mean cai, GC content and GenomeID
    genomeIDcol <- c(genomeIDcol, genomeID)
    meancaicol <- c(meancaicol, mean.cai)
    GCcontcol <- c(GCcontcol, GCcont)
    
  }
}

meancaiGCcontdf <- data.frame(genomeIDcol, meancaicol, GCcontcol,
                              stringsAsFactors = FALSE)
    
    
save(meancaiGCcontdf, file = "GCcontMeanCAI.RData")
load("GCcontMeanCAI.RData")



ggplot(data = meancaiGCcontdf, aes(x = meancaiGCcontdf[,2])) +
  geom_histogram(aes(y=..count../sum(..count..)), col = "darkgrey", fill = "black") +
  theme_bw(base_size = 13)+
  
  theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) +
  theme(legend.text = element_text(size = 17))  +
  #stat_bin(aes(y=..count../sum(..count..))) +
  #geom_histogram(aes(y=..count../sum(..count..)))  +
  labs(x = "mean CAI", 
       y = "number of genomes")


ggplot(data = meancaiGCcontdf, aes(x = meancaiGCcontdf[,3], y = meancaiGCcontdf[,2])) +
  geom_point(color = "blue", shape = 18) +
  theme_bw(base_size = 13)+
  theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) +
  theme(legend.text = element_text(size = 17))  +
  #stat_bin(aes(y=..count../sum(..count..))) +
  #geom_histogram(aes(y=..count../sum(..count..)))  +
  labs(x = "GC content (%)", 
       y = "mean CAI")



######______________________________For 1 genome______________________________######

genomeID <- "GCA_000003645"
cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
CDS.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
cai.data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)
seq.data <- read.csv(file = CDS.files, sep = ",", header = TRUE, as.is = TRUE)

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


mean.cai <- mean(cai.data[,2])


all.seq <- paste(as.matrix(seq.data)[,2], sep="", collapse="")
GC(all.seq)
seq.split <- strsplit(all.seq, "")[[1]]
GCcont <- GC(seq.split)*100
#nfreq <- table(seq.split)
#GCcont <- ((nfreq[2] + nfreq[3])/sum(nfreq))*100

