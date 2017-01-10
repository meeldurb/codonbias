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
w.data <- read.csv(file = "Reference_weight_tables_ENA/GCA_000003645.csv", 
                   header = FALSE, row.names = 1, as.is=TRUE) #as.is to keep the it as char
# only having numbers
w <- w.data[,1]

# retrieving domain data
domain.data <- read.csv(file = "Domain_data_ENA/GCA_000003645.csv", header = TRUE, 
                        as.is=TRUE) #as.is to keep the it as char

# renaming the columns
colnames(domain.data) <- c("domain_ID", "domain_beginpos", "domain_endpos", "CDS")
# checking how df looks like
str(domain.data)

####-----computing the postions of the domains with their corresponding postion in the CDS-----###

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
str(domain.data)

# iterate over all the CDS and "substring" the CDS that is associated to the protein domains
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
cai <- sapply(fasta.domains, cai, w = w)

# write the data to a file
write.table(cai, file = "CAI_GCA_000003645.csv", append = F, sep = "\t", row.names = names(cai), quote = F, col.names = F)

####################old script####################


# write this data to a file
write.table(CDS_dom_data, file = "20160428-Bpert_CDS_dom.csv", append = F, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#we converted the sequences to a fasta file with Python, to be able to be read by the CAI function
input_fasta <- read.fasta(file = "20160428-Bpert_CDS_dom.fasta")
w_data <- caitab$sc

# calculating the CAI for these sequences of XXX genome 
cai <- sapply(fasta.domains, cai, w = w_data)

# write to file 
# row.names are put with names(cai) to get the Pfam IDs in the file
write.table(cai, file = "20160428-cai_table_Bpert(GCA_001013565-1.LN849008).csv", append = F, sep = "\t", row.names = names(cai), quote = F, col.names = F)


