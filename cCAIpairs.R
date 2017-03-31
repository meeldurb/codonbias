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


genomeID <- "GCA_000003925"


# get files and open weight of codon pairs
w.files <- paste("codonpairs_weight/", genomeID, "_wcodpairs.csv", sep = "")
if (file.exists(w.files)){
  w.data <- read.csv(file = w.files, sep = ",", header = FALSE, as.is = TRUE)
  
}

# remove all codon pairs which start with stop codon
w.data <- w.data[order(w.data[,1]), ]
stop.pair.i <- grep("\\*_", w.data[,1])
length(w.data[,1]) - length(stop.pair.i)
w.pair <- w.data[-stop.pair.i, ]
length(w.pair[,1])


# make codon table for CAI calculations

# codon table 1 bacterial
aa1 <- getGeneticCode("SGC0")
# codon table 4 spiro/myco
aa4 <- getGeneticCode("SGC3")


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
codonpair.cai <- as.data.frame(aa.pairs)
codonpair.cai$codon.pairs <- codon.pairs

# aa4
codon.pairs <- NULL
aa.pairs <- NULL

for (first in names(aa4)){ 
  for (second in names(aa4)) {
    pcodons <- paste(first, second, sep="_")
    paa <- paste(aa4[first], aa4[second], sep="_")
    codon.pairs <- c(codon.pairs, pcodons)
    aa.pairs <- c(aa.pairs, paa)
  }
}

# making a dataframe from the codon pair table
aa4.codonpair.cai <- as.data.frame(aa.pairs)
aa4.codonpair.cai$codon.pairs <- codon.pairs



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



counter <- rep(0, length(codonpair.cai[,1]))
codonpair.cai$count <- counter

# counting RSCU of codon pairs
for (DNAseq in seqs2){
  # between every codon place an "_"
  sq <- gsub("(.{3})", "\\1 ", DNAseq)
  sq <- sub("\\s+$", "", sq)
  sq <- gsub(" ", "_", sq)
  # count the number of times a codon pair is present in the sequence
  for (codonpair in codonpair.cai[,2]){
    pair.count <- str_count(sq, pattern=codonpair)
    # when the count of codonpairs is above 0
    # it is first added to the df row of the associated codonpair
    # when there was already a count, the counts are added
    if (pair.count > 0){
      sel <- which(codonpair == codonpair.cai[,2])
      codonpair.cai[sel,3] <- sum(codonpair.cai[sel,3], pair.count)
    }
  }
}

codonpair.cai <- codonpair.cai[order(codonpair.cai[,1]), ]
countpairs.cai <- codonpair.cai[-stop.pair.i, ]

# change all NA's to 0 in cai df and w df
countpairs.cai[is.na(countpairs.cai)] <- 0
w.pair[is.na(w.pair)] <- 0


zero.threshold <- 0.000001
zero.to <- 0.0001
w.pair[w.pair[,3] < zero.threshold, 3] <- zero.to

# getting the sigma
sigma <- as.vector(countpairs.cai[,3]) %*% log(as.vector(w.pair[,3]))
exp(sigma/sum(w.pair[,3]))




####____________example calculating RSCU and weight_________####
seq1.1 <- tolower(strsplit(seq1, "")[[1]])
genomeID <- "GCA_000003645"
w.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
ordered.w <- w.data[with(w.data, order(V1)), ]
# only leaving numbers
w <- ordered.w[,2]
nncod1 <- uco(seq1.1)
sigma1 <- nncod1 %*% log(w)
sigma1
exp(sigma1/sum(nncod1))

