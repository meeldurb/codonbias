#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Compute weight table 
####  Purpose of script: The function to compute weight tables of 
####  data retrieved from SPARQL ENA database.
####  Only applicable for organisms with standard bacterial genetic code and 
####  Mycoplasma and Spiroplasma organism
#################################################################################




# # Packages that need to be installed
# source("http://bioconductor.org/biocLite.R")
# # ?BiocUpgrade
# biocLite("Biostrings")
# 
# # loading required libraries
# library("Biostrings")
# library("seqinr")
#



compute.weight <- function(seqs.df, genomeID){
  
  pasted.cdseqs <- paste(seqs.df, sep="", collapse="")
  # compute codon frequency
  codon.frequency <- trinucleotideFrequency(DNAString(pasted.cdseqs), as.prob=F, with.labels=T, as.array=F)
  
  
  # amino acid names of standard genetic code and myc/spiroplasma genetic code
  amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                   "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                   "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                   "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                   "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Stop", "Cys", "Trp", "Cys",
                   "Leu", "Phe", "Leu", "Phe")
  
  myc.spir.amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                            "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                            "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                            "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                            "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Trp", "Cys", "Trp", "Cys",
                            "Leu", "Phe", "Leu", "Phe")
  
  
  # Overwriting all the codons
  # changing codons that code for the same amino acid and dividing by the sum of it
  # Phenylalanine
  codon.frequency[c("TTT", "TTC")]<- codon.frequency[c("TTT", "TTC")]/sum(codon.frequency[c("TTT", "TTC")])
  #Leucine
  codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]<- codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]/sum(codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")])
  #Isoleucine
  codon.frequency[c("ATT", "ATC", "ATA")]<- codon.frequency[c("ATT", "ATC", "ATA")]/sum(codon.frequency[c("ATT", "ATC", "ATA")])
  #Valine
  codon.frequency[c("GTT", "GTC", "GTA", "GTG")]<- codon.frequency[c("GTT", "GTC", "GTA", "GTG")]/sum(codon.frequency[c("GTT", "GTC", "GTA", "GTG")])
  #Serine
  codon.frequency[c("TCT", "TCC", "TCA", "TCG")]<- codon.frequency[c("TCT", "TCC", "TCA", "TCG")]/sum(codon.frequency[c("TCT", "TCC", "TCA", "TCG")])
  #Proline
  codon.frequency[c("CCT", "CCC", "CCA", "CCG")]<- codon.frequency[c("CCT", "CCC", "CCA", "CCG")]/sum(codon.frequency[c("CCT", "CCC", "CCA", "CCG")])
  #threonine
  codon.frequency[c("ACT", "ACC", "ACA", "ACG")]<- codon.frequency[c("ACT", "ACC", "ACA", "ACG")]/sum(codon.frequency[c("ACT", "ACC", "ACA", "ACG")])
  #Alanine
  codon.frequency[c("GCT", "GCC", "GCA", "GCG")]<- codon.frequency[c("GCT", "GCC", "GCA", "GCG")]/sum(codon.frequency[c("GCT", "GCC", "GCA", "GCG")])
  #Tyrosine
  codon.frequency[c("TAT", "TAC")]<- codon.frequency[c("TAT", "TAC")]/sum(codon.frequency[c("TAT", "TAC")])
  #Histidine
  codon.frequency[c("CAT", "CAC")]<- codon.frequency[c("CAT", "CAC")]/sum(codon.frequency[c("CAT", "CAC")])
  #glutamine
  codon.frequency[c("CAA", "CAG")]<- codon.frequency[c("CAA", "CAG")]/sum(codon.frequency[c("CAA", "CAG")])
  #asparginine
  codon.frequency[c("AAT", "AAC")]<- codon.frequency[c("AAT", "AAC")]/sum(codon.frequency[c("AAT", "AAC")])
  #lysine
  codon.frequency[c("AAA", "AAG")]<- codon.frequency[c("AAA", "AAG")]/sum(codon.frequency[c("AAA", "AAG")])
  #aspartic acid
  codon.frequency[c("GAT", "GAC")]<- codon.frequency[c("GAT", "GAC")]/sum(codon.frequency[c("GAT", "GAC")])
  #glutamic acid
  codon.frequency[c("GAA", "GAG")]<- codon.frequency[c("GAA", "GAG")]/sum(codon.frequency[c("GAA", "GAG")])
  #cysteine
  codon.frequency[c("TGT", "TGC")]<- codon.frequency[c("TGT", "TGC")]/sum(codon.frequency[c("TGT", "TGC")])
  #arginine
  codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]<- codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]/sum(codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")])
  #serine
  codon.frequency[c("AGT", "AGC")]<- codon.frequency[c("AGT", "AGC")]/sum(codon.frequency[c("AGT", "AGC")])
  #glycine
  codon.frequency[c("GGT", "GGC", "GGA", "GGG")]<- codon.frequency[c("GGT", "GGC", "GGA", "GGG")]/sum(codon.frequency[c("GGT", "GGC", "GGA", "GGG")])
  #startcodon
  codon.frequency["ATG"]<- codon.frequency["ATG"]/codon.frequency["ATG"]
  
  #load all the genomeIDs and organisms for checking for Myco/Spiroplasma
  genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                   as.is=TRUE) #as.is to keep the it as char
  
  # Searching for Mycoplasma and Spiroplasma, which use other genetic codon table
  match.words <- c("Mycoplasma", "Spiroplasma")
  # i contains the indices where Myco/Spiro is found
  i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
  genomesMycoSpiro <- genome.and.organisms[i,1]
  # when the genome number of myc/spir is not found it will run the standard genetic code
  # else it will run the codon table of myc+spiroplasma
  if(!(genomeID %in% genomesMycoSpiro)) {
    #tryptophan (default)
    codon.frequency["TGG"]<- codon.frequency["TGG"]/codon.frequency["TGG"]
    #stopcodon (default)
    codon.frequency[c("TAA", "TAG", "TGA")]<- codon.frequency[c("TAA", "TAG", "TGA")]/sum(codon.frequency[c("TAA", "TAG", "TGA")])
    codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
    codon.table <- cbind(codon.frequency.table, amino.acids)
    sorted.codon.table <- codon.table[order(amino.acids),]
  } else { 
    #tryptophan (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
    codon.frequency[c("TGG", "TGA")]<- codon.frequency[c("TGG", "TGA")]/sum(codon.frequency[c("TGG", "TGA")])
    #stopcodon (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
    codon.frequency[c("TAA", "TAG")]<- codon.frequency[c("TAA", "TAG")]/sum(codon.frequency[c("TAA", "TAG")])
    codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
    codon.table <- cbind(codon.frequency.table, myc.spir.amino.acids)
    sorted.codon.table <- codon.table[order(myc.spir.amino.acids),]
  }
  
  # correct the codon frequencies for the maximum codon.frequency belonging to one amino acid
  tableFinalC=NULL; 
  tableFinalV=NULL; 
  tableFinalA=NULL;
  for (aa in unique(sorted.codon.table[,2])) { 
    #cat(aa, "\n")
    which.max <- which(sorted.codon.table[,2]==aa)
    tableFinalC <- c(tableFinalC, rownames(sorted.codon.table)[which.max])
    tableFinalV <- c(tableFinalV, sorted.codon.table[which.max,1]/max(sorted.codon.table[which.max,1] ))
    tableFinalA <- c(tableFinalA, rep(aa, length(which.max)))
  }
  final.codon.table <- data.frame(tableFinalC, tableFinalV, tableFinalA, stringsAsFactors = FALSE)
  colnames(final.codon.table) <- c("codon", "codon.frequency", "amino.acid")
  return(final.codon.table)
}

#genomeID <- "GCA_000003925"
#gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")

#gene.data <- read.csv(file = gene.files, header = TRUE, 
#                      as.is=TRUE) #as.is to keep the it as char

#w.table <- compute.weight(gene.data[,2], genomeID)

