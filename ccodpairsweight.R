#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Compute pairs weight table 
####  Purpose of script: computing codonpair weight tables function
####  Only applicable for organisms with standard bacterial genetic code and 
####  Mycoplasma and Spiroplasma organism
#################################################################################



# 
# # Packages that need to be installed
# source("http://bioconductor.org/biocLite.R")
# # ?BiocUpgrade
# biocLite("Biostrings")
# 
# # loading required libraries
# library("Biostrings")
# library("seqinr")
# 
# 
# genome.and.organisms <- read.csv(file = "test_genomes_ENA10.csv", header = FALSE, 
#                                  as.is=TRUE) #as.is to keep the it as char
# 
# outfolder <- "codonpairs_weight/"  
# if (!file.exists(outfolder))dir.create(outfolder)


    
compute.codpairs.weight <- function(w.seqs.df, genomeID){
  
#####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving codon tables ~~~~~~~~~~~~~~~~~~~~~~~~~#####
  
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
  aa1.codonpair.table <- as.data.frame(aa.pairs)
  aa1.codonpair.table$codon.pairs <- codon.pairs
  
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
  aa4.codonpair.table <- as.data.frame(aa.pairs)
  aa4.codonpair.table$codon.pairs <- codon.pairs
  
  
  
      # Searching for Mycoplasma and Spiroplasma, which use other genetic codon table
      match.words <- c("Mycoplasma", "Spiroplasma")
      # i contains the indices where Myco/Spiro is found
      i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
      genomesMycoSpiro <- genome.and.organisms[i,1]
      # when the genome number of myc/spir is not found it will run the standard genetic code
      # else it will run the codon table of myc+spiroplasma
      if(!(genomeID %in% genomesMycoSpiro)) {
        codonpair.table = aa1.codonpair.table
      } else { 
        codonpair.table = aa4.codonpair.table
      }
      # adding extra column to df for the counting of codonpairs
      counter <- rep(0, length(codonpair.table[,1]))
      codonpair.table$count <- counter
      # keeping seqcount to know at which DNAseq the loop is
      for (DNAseq in w.seqs.df){
        # between every codon place an "_"
        sq <- gsub("(.{3})", "\\1 ", DNAseq)
        sq <- sub("\\s+$", "", sq)
        sq <- gsub(" ", "_", sq)
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
          #print (codonpair.table[codonrow,2])
          freq <- count/totcount.aa
          codonpair.table[codonrow, 4] <- freq
          codonrow <- sum(codonrow, 1)
        }
      }
      # divide the frequencies by the maximum frequency of the aa pair
      final.cod <- NULL
      final.val <- NULL
      final.aa <- NULL
      for (aa in unique(codonpair.table[,1])){
        which.max <- which(codonpair.table[,1] == aa)
        final.cod <- c(final.cod, codonpair.table[which.max,2])
        final.val <- c(final.val, codonpair.table[which.max,4]/max(codonpair.table[which.max,4] ))
        final.aa <- c(final.aa, rep(aa, length(which.max)))
      }
      final.codonpair.table <- data.frame(final.aa, final.cod, final.val, stringsAsFactors = FALSE)
      
      return(final.codonpair.table)
    }

  

#genomeID <- "GCA_000003925"
#gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")

#gene.data <- read.csv(file = gene.files, header = TRUE, 
#                      as.is=TRUE) #as.is to keep the it as char

#w.codpairtable <- compute.codpairs.weight(genomeID)

  
