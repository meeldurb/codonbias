#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the CAI of the complete CDS of the genomes
##from the ENA namespace on blazegraph
###################################################################

# install the needed packages
install.packages("RCurl", repos="http://cran.rstudio.com/")
install.packages("XML", repos="http://cran.rstudio.com/")
# SPARQL needs RCurl and XML
install.packages("SPARQL", repos="http://cran.rstudio.com/")
# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")

# loading required library to execute SPARQL query
library("SPARQL")
# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")

outfolder <- "CAI_CDS/"  
if (!file.exists(outfolder))dir.create(outfolder)


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


# calculating the codon adaptation index for all genomes and its complete CDS

for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CAI_CDS.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    # Going through all genome numbers and retrieving their weight vectors and domain data
    w.files <- paste("Reference_weight_tables_ENA/", genomeID, ".csv", sep = "")
    if (file.exists(w.files)){
      w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
      # ordering on rownames, for it to be correctly used by cai function
      ordered.w <- w.data[with(w.data, order(V1)), ]
      # only leaving numbers
      w <- ordered.w[,2]
      gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
      gene.data <- read.csv(file = gene.files, header = TRUE, 
                              as.is=TRUE) #as.is to keep the it as char
      
     
      # writing to a file to be read by convert file
      write.table(gene.data, file = "tmp.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      # converting the sequences to fasta format with python file
      convertcmd <- "python Write_genecsvtofasta.py"
      system(convertcmd)
      
      # opening the written file and calculating the CAI
      fasta.genes <- read.fasta(file = "tmp.fasta")
      
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
      write.table(cai.df, file = fileout, append = F, sep = ",", quote = F, col.names = T, row.names = F)

      }
  }
}

#############_________________Function for computing cai from SPARQL data_________________#############

compute.cai <- function(seqs.df, genome_ID){

write.table(seqs.df, file = "tmpcai.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)

# converting the sequences to fasta format with python file
convertcmd <- "python Write_genecsvtofasta.py"
system(convertcmd)

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
