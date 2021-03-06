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



compute.cai <- function(seqs.df, genome_ID, w, tmpfile.seqs = "tmp.csv", tmpfile.fasta = "tmp.fasta"){
  
  write.table(seqs.df, file = tmpfile.seqs, append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # converting the sequences to fasta format with python file 
  # when it contains CDS the names of seqs consist of a link
  # when it contains domain seqs the names consist of pf_IDs
  
  #if (grepl("http://", seqs.df[1,1])){
  #convertcmd <- "python Write_genecsvtofasta.py"
  #}
  #if (grepl("pf", seqs.df[1,1])) {
  #convertcmd <- "python Write_csvtofasta.py"
  #}
  convertcmd <- paste("python write2fasta.py", tmpfile.seqs, tmpfile.fasta)
  system(convertcmd)
  # do something with ifelse statement when domains or when CDSs.
  
  # opening the written file and calculating the CAI
  fasta.genes <- read.fasta(file = tmpfile.fasta)
  
  
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
  rownames(cai.df) <- NULL
  return(cai.df)
}

# genomeID <- "GCA_000003645"
# gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
# gene.data <- read.csv(file = gene.files, header = TRUE, 
#                       as.is=TRUE) #as.is to keep the it as char
# seqs.df <- gene.data
# caivals <- compute.cai(gene.data, genomeID, w, "tmpcai.csv", "tmpcai.fasta")
