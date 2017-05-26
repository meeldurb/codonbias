#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for generating dataframe of all codons and genomes and 
## the relative adaptiveness
###################################################################


setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000008205"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char



####______________ create dataframe of all genomes with the codons and their relative adaptiveness ______________####
genomecount = 0
n = 0
c = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")  
  w.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "") 
  if (file.exists(w.files)){
    # Searching for Mycoplasma and Spiroplasma, which use other genetic codon table
    match.words <- c("Mycoplasma", "Spiroplasma")
    # i contains the indices where Myco/Spiro is found
    i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
    genomesMycoSpiro <- genome.and.organisms[i,1]
    # when the genome number of myc/spir is not found it will combine the genomeID in the dataframe
    # else the dataframe is too short (less stop codons are omitted)
    # so it will not complain anymore
    if(!(genomeID %in% genomesMycoSpiro)) {
      n <- n + 1
      w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
      # order on codon because of cai function
      #ordered.w <- w.data[with(w.data, order(w.data[,1])), ]
      # only leaving numbers
      w <- w.data[,2]
      if (genomecount == 0){
        codgendf <- data.frame(row.names=w.data[,1])
        codgendf <- cbind(codgendf, w)
        colnames(codgendf)[n] <- genomeID
        genomecount = genomecount + 1
      } else {
        codgendf <- cbind(codgendf, w)
        colnames(codgendf)[n] <- genomeID
      } 
    }else{
      print(paste("genomefile", "myco/spiro plasma"))
    }
  }else {
    print(paste("genomefile", genomeID, "does not exist"))
    c <- c + 1
  }
}
save(codgendf, file = "GenomeDataSet.RData")