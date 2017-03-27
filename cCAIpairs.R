#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing CAI of codon pairs 
####  Purpose of script: Function of computation of CAI of codon pairs
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
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/ccodpairsweight.R")



#############______________________ fake sequences ______________________###########
# counting the codon pairs
# fake sequence
seq1 <- "GGTGGTGGTGGCGGTGGCGGTGGCGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGG"
seq2 <- "GGAGGAGGAGGAGGGGGAGGGGGAGGGGGAGGAGGTGGAGGTGGAGGTGGAGGTGGAGGT"
seq3 <- "GGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGG"
seq4 <- "AACTGCAACTGCAACTGCAACTGCTACAGGTACAGGTACAGGTACAGGTACAGGTACAGG"
# how it should be
#seqs3 <- c(seq1, seq2)
# convert to df, because also read like this when loading the data
seqs <- as.data.frame(c(seq1, seq2, seq3, seq4))
seqs2 <- as.vector(seqs[,1])
#############______________________ fake sequences ______________________###########

zero.threshold <- 0.0001
zero.to <- 0.001

stop.pair.i <- grep("\\*_", w.codpairtable[,1])
length(w.codpairtable[,1]) - length(stop.pair.i)
w <- w.codpairtable[-stop.pair.i, ]
length(w[,1])
#w[which(w[,3] < zero.threshold), 3] <- zero.to
# after every letter put a space, then split on spaces 
# to get vector of chars
seq1 <- tolower(seq1)
sq <- gsub("(.{1})", "\\1 ", seq1)
sq1 <- strsplit(sq, " ")
sq1 <- sq1[[1]]

# calculate the RSCU of the sequence
nncod <- uco(sq1)
