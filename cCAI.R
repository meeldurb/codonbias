#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Compute weight table 
####  Purpose of script: The function to compute weight tables of 
####  data retrieved from SPARQL ENA database.
####  Only applicable for organisms with standard bacterial genetic code and 
####  Mycoplasma and Spiroplasma organism
#################################################################################




# Packages that need to be installed
# source("http://bioconductor.org/biocLite.R")
# ?BiocUpgrade
# biocLite("Biostrings")
# 
# # loading required libraries
# library("Biostrings")
# library("seqinr")



compute.cai <- function(seqs.df, genome_ID){
  
  write.table(seqs.df, file = "tmpcai.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # converting the sequences to fasta format with python file
  convertcmd <- "python Write_genecsvtofasta.py"
  system(convertcmd)
  # do something with ifelse statement when domains or when CDSs.
  
  # opening the written file and calculating the CAI
  fasta.genes <- read.fasta(file = "tmp.fasta")
  
  
  #load all the genomeIDs and organisms for checking for Myco/Spiroplasma
  genome.and.organisms <- read.csv(file = "genomes_ENA.csv", 
                                   header = FALSE, as.is=TRUE) 
  
  # Searching for Mycoplasma and Spiroplasma, using other weight tables
  match.words <- c("Mycoplasma", "Spiroplasma")
  # i contains the indices where Myco/Spiro is found
  i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
  genomesMycoSpiro <- genome.and.organisms[i,1]
  # when the genome number of myc/spir is not found it will run the default
  # else it will run with the codon table of myc+spiroplasma
  if(!(genomeID %in% genomesMycoSpiro)) {
    cai.output <- sapply(fasta.genes, cai, w = w, numcode = 1)
  } else { 
    cai.output <- sapply(fasta.genes, cai, w = w, numcode = 4)
  }
  # write the data to a file
  cai.df <- data.frame(cai.output)
  cai.df$gene_id <- names(cai.output)
  cai.df <- cai.df[,c(2,1)]
  return(cai.df)
}


compute.cai(gene.data, "GCA_000003645")