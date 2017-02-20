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
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")

# loading required libraries
library("Biostrings")
library("seqinr")

#loading my own functions
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000003645"



#weight table
w.files <- paste("Reference_weight_tables_ENA/", genomeID, ".csv", sep = "")
w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
# order on codon because of cai function
ordered.w <- w.data[with(w.data, order(w.data[,1])), ]
# only leaving numbers
w <- ordered.w[,2]

#CDS data
#cai.files <- paste("CAI_CDS/", genomeID, "_CAI_CDS.csv", sep = "")
#cai.data <- read.csv(file = cai.files, header = TRUE, as.is = TRUE)
gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
gene.data <- read.csv(file = gene.files, header = TRUE, 
                      as.is=TRUE) #as.is to keep the it as char


# compute CAI for the first round
cai.ini <- compute.cai(gene.data, genomeID, w, "tmpcai.csv", "tmpcai.fasta")

# sort on CAI value and take top 50
#sort.cai <- cai.data[order(-cai.data[,2]),]
ini.sort.cai <- cai.ini[order(-cai.ini[,2]),]
ini.top50 <- head(ini.sort.cai, 25)

# We retrieve only the CDS of the gene_IDs that are in the top 50
match.id <- as.vector(ini.top50[,1])
gene.match <- gene.data[gene.data[,1] %in% match.id, ] 

# then re-compute weight tables 
w.table <- compute.weight(gene.match[,2], genomeID)
ordered.w <- w.table[with(w.table, order(w.table[,1])), ]
# only leaving numbers
w <- ordered.w[,2]


# We re-calculate the CAI of all the gene_IDs 
cai.res <- compute.cai(gene.data, genomeID, w, "tmpcai.csv", "tmpcai.fasta")

# and take the top 50 again
res.sort.cai <- cai.res[order(-cai.res[,2]),] 
res.top50 <- head(res.sort.cai, 25)

#compare initial and result top 50
diff.count <- length(setdiff(ini.top50[,1], res.top50[,1]))

# if the difference between both lists is lower than 5, run the analysis again
# else save the resulting weight table
while (diff.count > 5){
  ini.top50 <- res.top50
  # We retrieve only the CDS of the gene_IDs that are in the top 50
  match.id <- as.vector(ini.top50[,1])
  gene.match <- gene.data[gene.data[,1] %in% match.id, ] 
  
  # then re-compute weight tables 
  w.table <- compute.weight(gene.match[,2], genomeID)
  ordered.w <- w.table[with(w.table, order(w.table[,1])), ]
  # only leaving numbers
  w <- ordered.w[,2]
  
  
  # We re-calculate the CAI of all the gene_IDs 
  cai.res <- compute.cai(gene.data, genomeID, w, "tmpcai.csv", "tmpcai.fasta")
  
  # and take the top 50 again
  res.sort.cai <- cai.res[order(-cai.res[,2]),] 
  res.top50 <- head(res.sort.cai, 25)
  
  #compare initial and result top 50
  diff.count <- length(setdiff(ini.top50[,1], res.top50[,1]))
  cat(paste("differences between tables is ", diff.count, "\n"))
  
} 

# creating the folder to save the data in
outfolder <- "Iterated_weight_tables_ENA/"  
if (!file.exists(outfolder))dir.create(outfolder)

fileout <- paste(outfolder, genomeID, "_it_weight25.csv", sep="")

write.table(w.table, file = fileout, append = FALSE, sep = ",", 
              row.names = FALSE, quote = FALSE, col.names = FALSE)







