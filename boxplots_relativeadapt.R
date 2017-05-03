#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script calculating and drawing boxplots from the 
##relative adaptiveness 
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
install.packages("plotly", repos="http://cran.rstudio.com/")
install.packages("ggplot2", repos= "http://cran.rstudio.com/")
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")


# loading required libraries

library("Biostrings")
library(ggplot2)
library(plotly)
library(ggplot2)
library("RColorBrewer")

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
# changing aa abbreviations from 1 to 3 letters
aatable[aatable == "A"] <- "Ala"
aatable[aatable == "R"] <- "Arg"
aatable[aatable == "N"] <- "Asn"
aatable[aatable == "D"] <- "Asp"
aatable[aatable == "C"] <- "Cys"
aatable[aatable == "E"] <- "Glu"
aatable[aatable == "Q"] <- "Gln"
aatable[aatable == "G"] <- "Gly"
aatable[aatable == "H"] <- "His"
aatable[aatable == "I"] <- "Ile"
aatable[aatable == "L"] <- "Leu"
aatable[aatable == "K"] <- "Lys"
aatable[aatable == "M"] <- "Met"
aatable[aatable == "F"] <- "Phe"
aatable[aatable == "P"] <- "Pro"
aatable[aatable == "S"] <- "Ser"
aatable[aatable == "T"] <- "Thr"
aatable[aatable == "W"] <- "Trp"
aatable[aatable == "Y"] <- "Tyr"
aatable[aatable == "V"] <- "Val"
aatable[aatable == "*"] <- "Stp"

# and order on aa and then on codon
ordered.aa <- aatable[with(aatable, order(aatable[,1], aatable[,2])), ]

# paste the aa names to the dataframe
codgendf <- data.frame(ordered.aa[,1], codgendf)
#codgendf <- data.frame(rownames(codgendf), codgendf)


mar.default <- c(4,3,1,1) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
par(mfrow = c(1,1))
par(cex.axis = 0.5)


# setting custom colors
palette(brewer.pal(length(unique(codgendf[,1])), "Accent"))(n)
boxplot(t(codgendf[,2:ncol(codgendf)]), col = codgendf[,1],
        xlab = "relative adaptiveness", 
        # main = "Boxplots of relative adaptiveness of 6000 genomes",
        horizontal = TRUE, las = 1)


### Draw plot
myplot <- ggplot(codgendf, 
                 aes(x = codgendf[,1], 
                     y = codgendf[,3:ncol(codgendf)], color=codgendf[,2]))+   
  #these commands create the plot, but nothing appears
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5))+
  theme_bw()+
  theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) +
  theme(legend.text = element_text(size = 10))  +
  #xlab(paste("PC1 (", format(pca.summary$importance[2,1]*100, digits = 2),"%)", sep = "")) +   
  #ylab(paste("PC2 (", format(pca.summary$importance[2,2]*100, digits = 2),"%)", sep = "")) +
  ggtitle(title) 

print(myplot)   #have a look at the plot 











