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





#####________________________________Making codonpair table________________________________#####
#codons <- words()
#aa <- translate(s2c(c2s(words())))
#names(aa) <- codons
#aa
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

# making a dataframe from the codon pair table
codonpair.table <- as.data.frame(aa.pairs)
codonpair.table$codon.pairs<- codon.pairs





#####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving ribosomal protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####

# open and read file with genomeIDs
genomeID <- "GCA_000003925"
genome.and.organisms <- read.csv(file = "test_genomes_ENA10.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# based on ribosomal protein domains
# We will go through all the genomes and the weight tables 
# will be written to a file
# creating a variable to store all the ribosomal proteins from all genomes

query.ribosomal.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?domain_id ?d_begin ?d_end ?CDS
WHERE {
VALUES ?sample { <http://gbol.life#xxx> }    
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:feature ?gene .
?gene gbol:transcript ?transcript .
?transcript gbol:sequence ?CDS .
?gene gbol:origin ?origin .
?origin <http://www.w3.org/ns/prov#provo:wasAttributedTo> ?annotation .
?annotation gbol:name "prodigal" .
?annotation gbol:version "2.6.3".
?transcript gbol:feature ?cds .
?cds gbol:protein ?protein .
?protein gbol:feature ?feature .
?feature gbol:provenance ?provenance .
?provenance gbol:signature_description ?domain_description .
?feature gbol:location ?location .
?location gbol:begin ?beginiri .
?beginiri gbol:position ?d_begin .
?location gbol:end ?endiri .
?endiri gbol:position ?d_end .
?feature gbol:origin <http://gbol.life#/interproscan/interproscan-5.21-60.0> .
?feature gbol:xref ?xref .
?xref gbol:PrimaryAccession ?domain_id .
?xref gbol:db <http://identifiers.org/pfam/> .
FILTER(regex(?domain_description, "ribosomal protein", "i"))
}
'

ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

outfolder <- "codonpairs_weight/"  
if (!file.exists(outfolder))dir.create(outfolder)

# looping through all the genomes
for (genomeID in genome.and.organisms[,1]) { 
  fileout <- paste(outfolder, genomeID, "_wcodpairs.csv", sep="")
  if (!file.exists(fileout)) {
    cat (genomeID, "\n")
    
    sub.query.ribosomal.seqs <- sub("xxx", genomeID, query.ribosomal.seqs)
    # running curl from command line
    curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.query.ribosomal.seqs,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
    curl <- gsub(pattern = "\n", replacement = " ", x = curl)
    system(curl)
    output.riboseqs <- rbind(sub.query.ribosomal.seqs, read.csv("tmp.txt", sep = "\t"))
    # slice off 1st row which contains NA values
    ribosomal.seqs.data <- output.riboseqs[-1,]
    #some genomes do not have ribosomal domain data, we should skip these
    if (length(ribosomal.seqs.data) == 0) {
      next 
      }
    # adding extra column to df for the counting of codonpairs
    counter <- rep(0, length(codonpair.table[,1]))
    codonpair.table$count <- counter
    # keeping seqcount to know at which DNAseq the loop is
    seqcount <- 1
    for (DNAseq in ribosomal.seqs.data[,4]){
      # between every codon place an "_"
      sq <- gsub("(.{3})", "\\1 ", DNAseq)
      sq <- sub("\\s+$", "", sq)
      sq <- gsub(" ", "_", sq)
      seqcount <- sum(seqcount, 1)
      print (seqcount)
      # count the number of times a codon pair is present in the sequence
      for (codonpair in codonpair.table[,2]){
        pair.count <- str_count(sq, pattern=codonpair)
        # when the count of codonpairs is above 0
        # it is first added to the df row of the associated codonpair
        # when there was already a count, the counts are added
        if (pair.count > 0){
          # take the codonpair from data that is the same as this one
          sel <- which(codonpair == codonpair.table[,2])
          codonpair.table[sel,3] <- sum(codonpair.table[sel,3], pair.count)
        }
      }
    }
    # add empty column to calculate frequencies
    codonpair.table$frequency <- rep(0, length(codonpair.table[,1]))
    # order the table on aa pairs, for later filling of frequencies
    codonpair.table <- codonpair.table[order(codonpair.table[,1]), ]
    # get the only one aa pair for looping
    unique.aa <- unique(codonpair.table[,1])
    # keep which codonrow the loop is at, to fill the frquencies
    codonrow <- 1
    # for every aa pair calculate the total of counts in every codon pair
    for (aa in unique.aa){
      aapair.count <- codonpair.table[which(codonpair.table[,1] == aa), ]
      totcount.aa <- sum(aapair.count[,3])
      print (aa)
      #print (totcount.aa)
      # then for every codonpair associated to an aa pair,
      # calculate the frequency of the codonpair (count/total count)
      for (count in aapair.count[,3]) {
        print (codonpair.table[codonrow,2])
        print (count)
        freq <- count/totcount.aa
        print (freq)
        codonpair.table[codonrow, 4] <- freq
        codonrow <- sum(codonrow, 1)
      }
    }
    write.table(codonpair.table, file = fileout, append = F, sep = ",", row.names = F, quote = F, col.names = F)
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

