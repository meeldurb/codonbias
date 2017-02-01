#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables 
####  Purpose of script: Re-computation of weight tables by selecting the top 50
####  domains with highest cai after each iteration is checked whether the top 50 
####  domains are comparable to previous iteration
#################################################################################

install.packages("fBasics", repos="http://cran.rstudio.com/")
# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")

# loading required libraries
library("Biostrings")
library("seqinr")
library("fBasics")
library("pwr")
library("effsize")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000003645"


#CDS data
domain.files <- paste("Domain_data_ENA/", genomeID, ".csv", sep = "")
domain.data <- read.csv(file = domain.files, header = TRUE, 
                        as.is=TRUE) #as.is to keep the it as char
colnames(domain.data) <- c("domain_ID", "domain_beginpos", "domain_endpos", "CDS")
keep <- c("domain_ID", "CDS")
domain.data.complete <- domain.data[keep]

#weight tables
w.files <- paste("Reference_weight_tables_ENA/", genomeID, ".csv", sep = "")
w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
ordered.w <- w.data[with(w.data, order(V1)), ]
# only leaving numbers
w <- ordered.w[,2]

# cai of all domains
cai.files <- paste("CAI_complete_ENA/", genomeID, "_CAI_complete.csv", sep = "")
cai.data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)
# sort on CAI value and take top 50
sort.cai <- cai.data[order(-V2),]
top50 <- head(sort.cai, 50)
# we only want to take the CDS of the exact same domains
# the domain ID should be the same but also starting position

match.pfid <- as.vector(top50[,1])
domain.match <- domain.data[domain.data$domain_ID %in% match.pfid,] 


#i <- grep(paste(match.pfid, collapse="|"), cai.data[,1])
#genomes.match <- cai.data[i,1]

match.pos <- as.vector(top50[,3])
pos.match <- domain.match[domain.match$domain_beginpos %in% match.pos,]

# Searching for Mycoplasma and Spiroplasma, using other weight tables
match.words <- c("Mycoplasma", "Spiroplasma")
# i contains the indices where Myco/Spiro is found
i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
genomesMycoSpiro <- genome.and.organisms[i,1]
# when the genome number of myc/spir is not found it will run the default
# else it will run the codon table of myc+spiroplasma
if(!(genomeID %in% genomesMycoSpiro)) {
  



# check which domains are duplicated
occurence <- rle(top50[,1])
duplicated.domains <- occurence$values[which(occurence$lengths!=1)]
for (dupl in duplicated.domains){
pos.dupl <- top50$V3[which(top50$V1 == dupl)]
print (pos.dupl)
}



for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CAI_complete.csv", sep="")
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
      keep <- c("domain_ID", "CDS")
      domain.data.complete <- domain.data[keep]
      
      # writing to a file to be read by convert file
      write.table(domain.data.complete, file = "tmp.csv", append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      