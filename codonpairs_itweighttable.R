#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables of codon pairs 
####  Purpose of script: Computation of iterative weight tables 
####  of codon pairs
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


setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open and read file with genomeIDs
#genomeID <- "GCA_000003925"
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


outfolder <- "codonpairs_weight/"  
if (!file.exists(outfolder))dir.create(outfolder)

# looping through all the genomes
for (genomeID in genome.and.organisms[,1]) { 
  fileout <- paste(outfolder, genomeID, "_wcodpairs.csv", sep="")
  if (!file.exists(fileout)) {
    cat (genomeID, "\n")
    w.codpairtable <- compute.codpairs.weight(genomeID)
    write.table(w.codpairtable, file = fileout, append = F, sep = ",", row.names = F, quote = F, col.names = F)
  }
}






#####___________________ example script ____________________####

#############______________________ fake sequences ______________________###########
# # counting the codon pairs
# # fake sequence
# seq1 <- "GGTGGTGGTGGCGGTGGCGGTGGCGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGG"
# seq2 <- "GGAGGAGGAGGAGGGGGAGGGGGAGGGGGAGGAGGTGGAGGTGGAGGTGGAGGTGGAGGT"
# seq3 <- "GGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGG"
# seq4 <- "AACTGCAACTGCAACTGCAACTGCTACAGGTACAGGTACAGGTACAGGTACAGGTACAGG"
# # how it should be
# #seqs3 <- c(seq1, seq2)
# # convert to df, because also read like this when loading the data
# seqs <- as.data.frame(c(seq1, seq2, seq3, seq4))
# seqs2 <- as.vector(seqs[,1])
#############______________________ fake sequences ______________________###########


# getting the files
genomeID <- "GCA_000003925"
gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
gene.data <- read.csv(file = gene.files, header = TRUE, 
                      as.is=TRUE) #as.is to keep the it as char


# adding extra column to df for the counting of codonpairs
counter <- rep(0, length(codonpair.table[,1]))
codonpair.table$count <- counter

seqcount <- 1
for (DNAseq in gene.data[,2]){
  # after every codon place an "_"
  sq <- gsub("(.{3})", "\\1 ", DNAseq)
  sq <- sub("\\s+$", "", sq)
  sq <- gsub(" ", "_", sq)
  seqcount <- sum(seqcount, 1)
  print (seqcount)
  for (codonpair in codonpair.table[,2]){
    pair.count <- str_count(sq, pattern=codonpair)
    if (pair.count > 0){
      # take the codonpair from data that is the same as this one
      sel <- which(codonpair == codonpair.table[,2])
      codonpair.table[sel,3] <- sum(codonpair.table[sel,3], pair.count)

    }
  }
}
codonpair.table$frequency <- rep(0, length(codonpair.table[,1]))
codonpair.table <- codonpair.table[order(codonpair.table[,1]), ]
unique.aa <- unique(codonpair.table[,1])
codonrow <- 1
for (aa in unique.aa){
  aapair.count <- codonpair.table[which(codonpair.table[,1] == aa), ]
  totcount.aa <- sum(aapair.count[,3])
  print (aa)
  #print (totcount.aa)
  for (count in aapair.count[,3]) {
    print (codonpair.table[meh,2])
    print (count)
    freq <- count/totcount.aa
    print (freq)
    codonpair.table[codonrow, 4] <- freq
    meh <- sum(codonrow, 1)
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

