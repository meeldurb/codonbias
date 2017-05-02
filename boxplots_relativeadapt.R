#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script calculating and drawing boxplots from the 
##relative adaptiveness 
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
install.packages("plotly", repos="http://cran.rstudio.com/")
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")


# loading required libraries

library("Biostrings")
library(ggplot2)
library(plotly)

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000008185"
genomeID <- "GCA_000008205"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genome.and.organisms <- genome.and.organisms[,1:2]


####________________________ Retrieve data of all the codon weights ________________________####

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

load("GenomeDataSet.RData")


####_______ separating the codons based on aa  _______####

aa1 <- getGeneticCode("SGC0")

unique.aa <- unique(codpairs[,1])

# keep which plot it is drawing, to keep all plots
plot = 1
# for every aa pair calculate the total of counts in every codon pair
for (aa in unique.aa){
  aapair.count <- codpairs[which(codpairs[,1] == aa), ]
  print (aa)
  codonrows <- NULL
  for (codpair in aapair.count[,2]){
    codpair.count <- which(rownames(codgendf) == codpair )
    codonrows <- c(codonrows, codpair.count)
  }
  filename <- paste(outfolder, "boxplot_codpairs_", plot, ".jpg", sep = "")
  if (plot == 1){
    jpeg(file=filename)
    # setting the format of screen, so that codonpairs are not drawn out of screen
    mar.default <- c(2,4,2,2) + 0.1
    par(mar = mar.default + c(0, 2, 0, 0))
    par(mfrow = c(3,2)) 
  }
  boxplot(t(codgendf[codonrows, ]), col = rainbow(20),
          xlab = "relative adaptiveness", 
          main = paste("Boxplots of codon pairs of aa pair ", aa, sep = "" ),
          horizontal = T, las = 1)
  plot = sum(plot, 1)
  if (plot %% 6 == 0 ){
    dev.off()
    jpeg(file=filename)
    # setting the format of screen, so that codonpairs are not drawn out of screen
    mar.default <- c(2,4,2,2) + 0.1
    par(mar = mar.default + c(0, 2, 0, 0))
    par(mfrow = c(3,2))
  }
}


