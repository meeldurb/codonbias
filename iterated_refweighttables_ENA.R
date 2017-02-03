#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables 
####  Purpose of script: Re-computation of weight tables by selecting the top 50
####  domains with highest cai after each iteration is checked whether the top 50 
####  domains are comparable to previous iteration
#################################################################################

# install needed packages
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")

# loading required libraries
library("Biostrings")
library("seqinr")

#loading my own functions
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")




setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000003645"


#CDS data
gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
gene.data <- read.csv(file = gene.files, header = TRUE, 
                        as.is=TRUE) #as.is to keep the it as char
colnames(gene.data) <- c("gene_ID", "CDS")

#weight table
w.files <- paste("Reference_weight_tables_ENA/", genomeID, ".csv", sep = "")
w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
# order on codon because of cai function
ordered.w <- w.data[with(w.data, order(V1)), ]
# only leaving numbers
w <- ordered.w[,2]

#cai values of the first round
cai.files <- paste("CAI_CDS/", genomeID, "_CAI_CDS.csv", sep = "")
cai.data <- read.csv(file = cai.files, header = TRUE, as.is = TRUE)





# sort on CAI value and take top 50
sort.cai <- cai.data[order(-cai.data$cai.output),]
ini.top50 <- head(sort.cai, 50)

# We retrieve only the CDS of the gene_IDs that are in the top 50
match.id <- as.vector(ini.top50[,1])
gene.match <- gene.data[gene.data$gene_ID %in% match.id,] 

codon.table <- compute.weight(gene.match[,2], "GCA_000003645")


# We re-calculate the CAI of all the gene_IDs 



write.table(gene.data, file = "tmpcai.csv", append = FALSE, sep = ",", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)




