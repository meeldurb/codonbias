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

genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# open the files with biased and unbiased genomes for colouring the plot
biased.genomes <- read.csv(file = "biased.genomes.csv", header = FALSE, 
                           as.is=TRUE)
biased.genomes <- as.vector(as.matrix(biased.genomes))

unbiased.genomes <- read.csv(file = "unbiased.genomes.csv", header = FALSE, 
                             as.is=TRUE)
unbiased.genomes <- as.vector(as.matrix(unbiased.genomes))

genomecount = 0
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
    xlim = c(10, 90)
    ylim = c(0.1, 1)
    # then plot the GC content against the mean CAI
      if(genomeID %in% biased.genomes){
        print ("biased genome")
        col = "blue"
        } else if (genomeID %in% unbiased.genomes){
          print ("unbiased genome")
        col = "red"
        } else {
        col = "grey"
        }
    if (genomecount == 0){    
      plot(GCcont, mean.cai, type = "p", xlim = xlim, ylim = ylim,
           pch = 18, col = col,
           main = "Average CAI vs. GC content",
           xlab = "GC content (%)", 
           ylab = "Mean CAI")
      grid(NULL, NULL, lty = 6, col = "cornsilk2")
      genomecount = genomecount + 1
      } else {
      points(GCcont, mean.cai, pch = 18, col = col)
      }
    legend("topright" ,c("biased", "unbiased", "unknown"), cex=1.5, pch=18,
           col=c("blue", "red", "grey") , bty="n")
  }
}



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

