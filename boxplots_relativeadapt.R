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

# get the codontable 1 and convert it to a dataframe
aa1 <- getGeneticCode("SGC0")
aatable <- data.frame(keyName = aa1, value = names(aa1), 
                      row.names = NULL, stringsAsFactors = FALSE)
colnames(aatable) <- c("aa", "codon")
aatable[aatable == "*"] <- "Stp"
order.aa <- aatable[with(aatable, order(aatable[,1], aatable[,2])), ]

# get the positons of the codons belonging to which amino acid
# in the big dataframe with all the genomes

ordered.aa <- NULL
for (codon in rownames(codgendf)){
  print (codon)
  position <- which(codon == aatable[,2])
  ordered.aa <- c(ordered.aa, aatable[position,1])
  print(aatable[position,1])
  
}

mar.default <- c(4,2,2,1) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
par(mfrow = c(1,1))
par(cex.axis = 0.5)

boxplot(t(codgendf), col = rainbow(20),
        xlab = "relative adaptiveness", 
        main = "Boxplots of relative adaptiveness of 6000 genomes",
        horizontal = TRUE, las = 1)

for (aa in unique.aa){
  # is drawing boxplots in one frame per amino acid. Want it in one screen
  aa.count <- aatable[which(aatable[,1] == aa), ]
  print (aa)
  cod.rows <- NULL
  for (codon in aa.count[,2]){
    cod.count <- which(rownames(codgendf) == codon )
    cod.rows <- c(cod.count,cod.rows)
  }
  boxplot(t(codgendf[cod.rows, ]), col = rainbow(20),
          xlab = "relative adaptiveness", 
          main = paste("Boxplots of codon weights", aa,
          horizontal = T, las = 1)
  )
}


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


