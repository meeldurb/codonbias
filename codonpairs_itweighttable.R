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
# for strings
install.packages("stringr", repos="http://cran.rstudio.com/")

# loading required libraries
library("Biostrings")
library("seqinr")
library("stringr")

#loading my own functions
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")



setwd("~/Documents/Master_Thesis_SSB/git_scripts")


codons <- words()
aa <- translate(s2c(c2s(words())))
names(aa) <- codons
aa
# codon table 1 bacterial
aa1 <- getGeneticCode("SGC0")
# codon table 4 spiro/myco
aa4 <- getGeneticCode("SGC3")



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

codonpair.table <- aa.pairs
names(codonpair.table) <- codon.pairs
codonpair.table


# fake sequence
seq1<- "AACTGCGATGATACAGGGATACCGATAGACTAGATCCAGGTAACGATAAAGCGAGAGGCC"
seq2 <- "AACTGCAATGATACAGGGATACGGGCTAACTAGATCCAGGTAACGATAAACATTGGAGGC"
seqs <- c(seq1, seq2)
for (DNAseq in seqs){
  sq <- gsub("(.{3})", "\\1 ", DNAseq)
  sq <- sub("\\s+$", "", sq)
  sq <- gsub(" ", "_", sq)
  print(sq)
  for (codonpair in names(codonpair.table)){
    #print (codonpair)
    meh <- str_count(sq, pattern=codonpair)
    if (meh > 0){
      print (codonpair)
      print (meh)
    }
  }
}

# counting occurence in string
#str_count(sq, "AAC_TGC")




#fakedata
data <- sample(1:100, 25)
names(data) <- pairsCodons
data
AAp <- unique(tranlPairs)

for (aap in AAp){
  sel <- names(which(tranlPairs==aap))
  data[sel]
}

AAp
pairsCodons
data

