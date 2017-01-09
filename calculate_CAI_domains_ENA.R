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



