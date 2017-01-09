##Author: Melanie van den Bosch
##Script for calculating the CAI of all the genomes and their domains
##from the ENA namespace on blazegraph


# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")
# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr")

# loading required library to execute SPARQL query
library("SPARQL")
# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# Retrieving weight vectors
w <- read.csv(file = "Reference_weight_tables_ENA/GCA_000003645.csv", header = FALSE, 
              as.is=TRUE) #as.is to keep the it as char


# retrieving domain data
domain.data <- read.csv(file = "Domain_data_ENA/GCA_000003645.csv", header = TRUE, 
                        as.is=TRUE) #as.is to keep the it as char
# renaming the columns
colnames(domain.data) <- c("domain_ID", "domain_beginpos", "domain_endpos", "CDS")
#checking how df looks like
str(domain.data)

####-----computing the postions of the domains with their corresponding postion in the CDS-----###

for (row in domain.data) 
{
  CDS_begin <- 3*(domain.data[2]-1)+1
  CDS_end <- 3*(domain.data[3])
}



# renaming the columns
colnames(CDS_begin) <- NULL
colnames(CDS_end) <- NULL

# binding these columns to already existing dataframe
domain.data$CDS_domain_beginpos <- as.integer(CDS_begin[,1])
domain.data$CDS_domain_end <- as.integer(CDS_end[,1])
str(domain.data)

# iterate over all the CDS and "substring" the CDS that is associated to the protein domains
for (CDS in domain.data[,4]){
  dom_seq <- substr(CDS, domain.data[domain.data$CDS == CDS, 5], 
                    domain.data[domain.data$CDS == CDS, 6])
  domain.data$CDS_domain <- dom_seq # for now it is incorrect, it adds only one sequence multiple times
}

####################old script####################

# creating a dataframe with the Pfam_ID, CDS_begin, CDS_end and the CDS_seq
CDS_dom <- cbind(Pfam_ID, CDS_begin_intTBL, CDS_end_intTBL, CDS_seq)

# setting the variables to NULL
seq_subTBL <- NULL
seq_sub <- NULL
seq <- NULL

# iterate over all the sequences in the dataframe and substracting the CDS that is associated to the protein domains
for (seq in CDS_dom[,4]){
  
  seq_sub <- substr(seq, CDS_dom[CDS_dom$CDS_seq == seq, 2], 
                    CDS_dom[CDS_dom$CDS_seq == seq, 3])
  #print (seq_sub)
  seq_subTBL <- append(seq_subTBL, seq_sub)
  seq_sub <- NULL
  
}

# Can be replaced, not working yet.
# seq_pos <- apply(CDS_dom[,4], 1, substr, CDS_dom[,2], CDS_dom[,3])

seq_subTBL <- as.data.frame(x=seq_subTBL)
colnames(seq_subTBL) <- "Truncated_Seqs"

# The dataframe now contains the Pfam ID, the start and ending position of the domain in the gene sequence
# the original CDS and the calculated CDS which are the sequences that are associated to the protein domains.
CDS_dom_data <- cbind(CDS_dom, seq_subTBL)

# write this data to a file
write.table(CDS_dom_data, file = "20160428-Bpert_CDS_dom.csv", append = F, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#we converted the sequences to a fasta file with Python, to be able to be read by the CAI function
input_fasta <- read.fasta(file = "20160428-Bpert_CDS_dom.fasta")

# calculating the CAI for these sequences of XXX genome 
cai <- sapply(input_fasta, cai, w = w_data)

# write to file 
# row.names are put with names(cai) to get the Pfam IDs in the file
write.table(cai, file = "20160428-cai_table_Bpert(GCA_001013565-1.LN849008).csv", append = F, sep = "\t", row.names = names(cai), quote = F, col.names = F)


