#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the CAI of all the genomes and their domains
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

outfolder <- "CAI_domains_ENA/"  
if (!file.exists(outfolder))dir.create(outfolder)


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


# calculating the codon adaptation index for all genomes and its domains

for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CAI.csv", sep="")
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
      domain.files <- paste("Domain_data_ENA/", genomeID, ".csv", sep = "")
      domain.data <- read.csv(file = domain.files, header = TRUE, 
                              as.is=TRUE) #as.is to keep the it as char
      # renaming the columns
      colnames(domain.data) <- c("domain_ID", "domain_beginpos", "domain_endpos", "CDS")
      
      # computing the postions of the domains with their corresponding postion in the CDS
      for (row in domain.data) 
      {
        CDS_begin <- 3*(domain.data[2]-1)+1
        CDS_end <- 3*(domain.data[3])
      }
      
      # renaming the columns
      colnames(CDS_begin) <- ""
      colnames(CDS_end) <- ""
      
      # binding these columns to already existing dataframe
      # and converting to integers for later calculation
      domain.data$CDS_domain_beginpos <- as.integer(CDS_begin[,1])
      domain.data$CDS_domain_end <- as.integer(CDS_end[,1])
      
      # Substracting the CDS part that is associated to the protein domain
      # from columns 5 and 6 the positions of the begin and end of the 
      # CDS associated to protein domain are derived
      # when CDS are equal, take the integer from that column
      for (CDS in domain.data[4]){
        dom_seq <- substr(CDS, domain.data[domain.data$CDS == CDS, 5], 
                          domain.data[domain.data$CDS == CDS, 6])
        
      }
      domain.data$CDS_domain <- dom_seq 
      
      # shrinking the data to reduce computational time
      keep <- c("domain_ID", "CDS_domain")
      domain.data.CDS <- domain.data[keep]
      
      write.table(domain.data.CDS, file = "tmp.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      # converting the sequences to fasta format with python file
      convertcmd <- "python Write_csvtofasta.py"
      system(convertcmd)
      
      # opening the written file and calculating the CAI
      fasta.domains <- read.fasta(file = "tmp.fasta")
      
      # Searching for Mycoplasma and Spiroplasma, using other weight tables
      match.words <- c("Mycoplasma", "Spiroplasma")
      # i contains the indices where Myco/Spiro is found
      i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
      genomesMycoSpiro <- genome.and.organisms[i,1]
      # when the genome number of myc/spir is not found it will run the default
      # else it will run with the codon table of myc+spiroplasma
      if(!(genomeID %in% genomesMycoSpiro)) {
        
      cai.output <- sapply(fasta.domains, cai, w = w, numcode = 1)
      } else { 
        cai.output <- sapply(fasta.domains, cai, w = w, numcode = 4)
      }
      # write the data to a file
      write.table(cai.output, file = fileout, append = F, sep = ",", row.names = names(cai.output), quote = F, col.names = F)

      }

  
}
}


####################_____________ calculating cai for single file _______________________################

# # Retrieving weight vectors
# w.data <- read.csv(file = "Reference_weight_tables_ENA/GCA_000003645.csv", 
#                    header = FALSE, row.names = 1, as.is=TRUE) #as.is to keep the it as char
# 
# # ordering on rownumbers for it to look the same as caitab
# w.data$names <-rownames(w.data)
# ordered.w <- w.data[with(w.data, order(names)), ]
# 
# # only having numbers
# w <- w.data[,1]
# w.ordered <- ordered.w[,1]
# 
# # retrieving domain data
# domain.data <- read.csv(file = "Domain_data_ENA/GCA_000003645.csv", header = TRUE, 
#                         as.is=TRUE) #as.is to keep the it as char
# 
# # renaming the columns
# colnames(domain.data) <- c("domain_ID", "domain_beginpos", "domain_endpos", "CDS")
# # checking how df looks like
# str(domain.data)
# 
# ####-----computing the postions of the domains with their corresponding postion in the CDS-----###
# 
# for (row in domain.data) 
# {
#   CDS_begin <- 3*(domain.data[2]-1)+1
#   CDS_end <- 3*(domain.data[3])
# }
# 
# # renaming the columns
# colnames(CDS_begin) <- ""
# colnames(CDS_end) <- ""
# 
# # binding these columns to already existing dataframe
# # and converting to integers for later calculation
# domain.data$CDS_domain_beginpos <- as.integer(CDS_begin[,1])
# domain.data$CDS_domain_end <- as.integer(CDS_end[,1])
# str(domain.data)
# 
# # iterate over all the CDS and "substring" the CDS that is associated to the protein domains
# for (CDS in domain.data[4]){
#   dom_seq <- substr(CDS, domain.data[domain.data$CDS == CDS, 5], 
#                     domain.data[domain.data$CDS == CDS, 6])
#   
# }
# domain.data$CDS_domain <- dom_seq 
# 
# # shrinking the data to reduce computational time
# keep <- c("domain_ID", "CDS_domain")
# domain.data.CDS <- domain.data[keep]
# 
# write.table(domain.data.CDS, file = "tmp.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
# 
# # converting the sequences to fasta format with python file
# convertcmd <- "python Write_csvtofasta.py"
# system(convertcmd)
# 
# # opening the written file and calculating the CAI
# fasta.domains <- read.fasta(file = "tmp.fasta")
# cai.output <- sapply(fasta.domains, cai, w = w)
# cai.output.ord <- sapply(fasta.domains, cai, w = w.ordered, numcode = 1)
# 
# # write the data to a file
# write.table(cai.output, file = "CAI_GCA_000003645.csv", append = F, sep = "\t", row.names = names(cai.output), quote = F, col.names = F)
# write.table(cai.output.ord, file = "CAI_GCA_000003645ord2.csv", append = F, sep = "\t", row.names = names(cai.output.ord), quote = F, col.names = F)
# 
# # gives different output.
# # found that cai function finds the indexpositions of codons that should be excluded from the analysis based on alphabetical order
