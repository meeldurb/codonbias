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

outfolder <- "CAI_domains_ENA/"  
if (!file.exists(outfolder))dir.create(outfolder)


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# Retrieving weight vectors
w.data <- read.csv(file = "Reference_weight_tables_ENA/xxx.csv", 
                   header = FALSE, row.names = 1, as.is=TRUE) #as.is to keep the it as char


# retrieving domain data
domain.data <- read.csv(file = "Domain_data_ENA/xxx.csv", header = TRUE, 
                        as.is=TRUE) #as.is to keep the it as char

for (genomeID in genome.and.organisms[,1]) { 
  fileout <- paste(outfolder, genomeID, "_CAI.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    sub.w.data <- sub("xxx", genomeID, w.data)
    print (sub.w.data)
    w <- sub.w.data[,1]
    sub.domain.data <- sub("xxx", genomeID, domain.data)
}
}


####################_____________ calculating cai for single file _______________________################

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

