#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the CAI of all the genomes and their domains
##from the ENA namespace on blazegraph
###################################################################


# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")
# to use string methods
install.packages("stringr", repos="http://cran.rstudio.com/")

# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")
# loading required library to use string methods
library("stringr")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")

outfolder <- "CAI_inter_domains_ENA/"  
if (!file.exists(outfolder))dir.create(outfolder)


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "test_genomes_ENA10.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


# calculating the codon adaptation index for all genomes and its domains

for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CAI_inter.csv", sep="")
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
      keep <- c("domain_ID", "CDS", "CDS_domain")
      domain.data.CDS <- domain.data[keep]
      
      for (i in 1:length(domain.data.CDS[1:5,1])){
        inter.dom.CDS <- gsub(pattern = test[i,3], replacement="", x=test[i,2])
        str(inter.dom.CDS)
        }
      }
  }
}




domain.data.CDS$CDS_inter_domain <- inter.dom.CDS
x <- "Do you wish you were Mr. Jones?"
y <- strsplit(x, ". ")
z <- paste(y, sep="") 

txt <- c("arm","foot","lefroo", "bafoobar")
if(length(i <- grep("foo", txt, invert = TRUE)))
  cat("'foo' appears at least once in\n\t", txt, "\n")
i # 2 and 4
pasted <- paste(txt[i], collapse='')



