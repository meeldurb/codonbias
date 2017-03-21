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





#####________________________________Example script________________________________#####
codons <- words()
aa <- translate(s2c(c2s(words())))
names(aa) <- codons
aa
# codon table 1 bacterial
aa1 <- getGeneticCode("SGC0")
# codon table 4 spiro/myco
aa4 <- getGeneticCode("SGC3")


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

codonpair.table <- as.data.frame(aa.pairs)
codonpair.table$codon.pairs<- codon.pairs
#length(codonpair.table[,2])


# counting the codon pairs
# fake sequence
seq1 <- "TACAGGTACAGGTACAGGTACAGGTACAGGTACAGGTACAGGTACAGGTACAGGTACAGG"
seq2 <- "AAGTGCAATGATACAGGGATACGGGCTAACTAGATCCAGGTAACGATAAACATTGGAGGC"
seq3 <- "AACTGCAACTGCAACTGCAACTGCAACTGCAACTGCAACTGCAACTGCAACTGCAACTGC"
# how it should be
#seqs3 <- c(seq1, seq2)
# convert to df, because also read like this when loading the data
seqs <- as.data.frame(c(seq1, seq2, seq3))
seqs2 <- as.vector(seqs[,1])

data <- rep(0, length(codonpair.table[,1]))
codonpair.table$count <- data

for (DNAseq in seqs2){
  sq <- gsub("(.{3})", "\\1 ", DNAseq)
  sq <- sub("\\s+$", "", sq)
  sq <- gsub(" ", "_", sq)
  print(sq)
  for (codonpair in codonpair.table[,2]){
    #print (codonpair)
    pair.count <- str_count(sq, pattern=codonpair)
    if (pair.count > 0){
      # take the codonpair from data that is the same as this one
      sel <- which(codonpair == codonpair.table[,2])
      codonpair.table[sel,3] <- sum(codonpair.table[sel,3], pair.count)
    }
  }
}

# check whether codons are saved in dataframe
findcod <- which("AAC_TGC" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("TGC_AAC" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("TAC_AGG" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("AGG_TAC" == codonpair.table[,2])
codonpair.table[findcod,3]






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

