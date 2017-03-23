#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the mean CAI vs the 
##GC content of all genomes grouped on phyla that are known
##which polIII isoform they contain
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



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

#genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# the 3 different groups of polIII isoforms and its phyla belonging to it
org.dnaE1 <- c("Actinobacteria", "Aquificae", "Bacteroidetes", "Chlorobi", "Chlamydiae", "Verrucomicrobia",
           "Chloroflexi", "Cyanobacteria", "Deinococcus-Thermus", "Proteobacteria", "Spirochaetes")

org.dnaE2.dnaE1 <- c("Fibrobacteres", "Acidobacteria", "Planctomycetes", "Proteobacteria")

org.polC.dnaE3 <- c("Firmicutes", "Fusobacteria", "Thermotogae")


# Taking the information on phyla from the gold db
gold.data= read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)

# taking the genomeIDs that have phyla information
phyla <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data)=="NCBI.Phylum")]))
phyla <- na.omit(phyla)
length(phyla$NCBI.Phylum)
# check which groups are there and how large they are
unique(phyla$NCBI.Phylum)
table(phyla[,2])

# making a list of which 
group.dnaE1 <- phyla[which(phyla$NCBI.Phylum %in% org.dnaE1), 1]
length(group.dnaE1)
group.dnaE2.dnaE1 <- phyla[which(phyla$NCBI.Phylum %in% org.dnaE2.dnaE1), 1]
length(group.dnaE2.dnaE1)
group.polC.dnaE3 <- phyla[which(phyla$NCBI.Phylum %in% org.polC.dnaE3), 1]
length(group.polC.dnaE3)

length(group.dnaE1) + length(group.dnaE2.dnaE1) + length(group.polC.dnaE3)

genomecount = 0
pdf("GCvsCAI_plot_itcount.pdf")
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
    xlim = c(20, 80)
    ylim = c(0.35, 0.8)
    if (genomeID %in% group.dnaE1){
      col = "blue"
    } else if (genomeID %in% group.dnaE2.dnaE1){
        col = "red"
    } else if(genomeID %in% group.polC.dnaE3) {
      col = "green"
    } else {
      col = "grey"
    }
    # then plot the GC content against the mean CAI
      # #if(genomeID %in% biased.genomes){
      #   print ("biased genome")
      #   col = "blue"
      #   } else if (genomeID %in% unbiased.genomes){
      #     print ("unbiased genome")
      #   col = "red"
      #   } else {
      #   col = "grey"
      #   }
    if (genomecount == 0){  
      plot(GCcont, mean.cai, col=col, 
                   type = "p", xlim = xlim, ylim = ylim,
           pch = 18, main = "Average CAI vs. GC content",
           xlab = "GC content (%)", 
           ylab = "Mean CAI")
      grid(NULL, NULL, lty = 6, col = "cornsilk2")
      legend("bottomright" ,c("dnaE1", "dnaE2/dnaE1", "polC/dnaE3", "Other"), cex=1.5, pch=18,
            col=c("blue", "red", "green", "grey") , bty="n")
      genomecount = genomecount + 1
      } else {
      points(GCcont, mean.cai, col= col, pch = 18)
      }
  }
}

dev.off()



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

