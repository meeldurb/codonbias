#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables 
####  Purpose of script: Re-computation of weight tables by selecting the top 25
####  domains with highest cai after each iteration is checked whether the top 25 
####  domains are comparable to previous iteration
#################################################################################

# install needed packages
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")

# loading required libraries
library("Biostrings")
library("seqinr")

#loading my own functions
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")



setwd("~/Documents/Master_Thesis_SSB/git_scripts")


codons <- words()
aa <- translate(s2c(c2s(words())))
names(aa) <- codons
aa
# codon table 1 bacterial
aa1 <- getGeneticCode("SGC1")
# codon table 4 spiro/myco
aa4 <- getGeneticCode("SGC4")




pairsCodons <- NULL
pairsAA <- NULL

for ( first in names(transl)){ 
  for (second in names(transl)) {
    pcodons <- paste(first, second, sep="_")
    paa <- paste(transl[first], transl[second], sep="_")
    pairsCodons <- c( pairsCodons,pcodons)
    pairsAA <- c(pairsAA, paa)
  }
}
pairsCodons
pairsAA
tranlPairs <- pairsAA
names(tranlPairs) <- pairsCodons

#fakedata
data <- sample(1:100, 25)
names(data) <- pairsCodons
AAp <- unique(tranlPairs)

for (aap in AAp){
  sel <- names(which(tranlPairs==aap))
  data[sel]
}

AAp
pairsCodons
data

