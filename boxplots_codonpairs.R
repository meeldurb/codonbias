#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script calculating and drawing boxplots from the codonpair 
##weight tables of all the genomes
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
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


####________________________ Retrieve data of all the codon pair weights ________________________####

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

load("codonpairsGenomeDataSet.RData")
#codgendf[is.na(codgendf)] <- 0

#####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving codon tables ~~~~~~~~~~~~~~~~~~~~~~~~~#####
# codon table 1 bacterial
aa1 <- getGeneticCode("SGC0")

# form the codonpair table for codontable 1 and 4 seperately
# aa1
# making the codonpair table with codons and aa
codon.pairs <- NULL
aa.pairs <- NULL


for (first in names(aa1)){ 
  for (second in names(aa1)) {
    pcodons <- paste(first, second, sep="_")
    paa <- paste(aa1[first], aa1[second], sep="_")
    codon.pairs <- c(codon.pairs, pcodons)
    aa.pairs <- c(aa.pairs, paa)
  }
}

# making a dataframe from the codon pair table
codonpair.table <- as.data.frame(aa.pairs)
codonpair.table$codon.pairs <- codon.pairs

# order on codon to get better readability
codonpair.table <- codonpair.table[with(codonpair.table, order(codonpair.table[,1])), ]
# removing all codon pairs which start with stop codon
# for weight tables
stop.pair.i <- grep("\\*_", codonpair.table[,1])
length(codonpair.table[,1]) - length(stop.pair.i)
codpairs <- codonpair.table[-stop.pair.i, ]
length(codpairs[,1])


####_______ splitting the codonpairs based on aa pairs per boxplot _______####

outfolder <- "codonpairs_boxplots/"  
if (!file.exists(outfolder))dir.create(outfolder)

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


