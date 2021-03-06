#############################################################################
####  Author: Melanie van den Bosch
####  Script for retrieving all genetic sequences of all genomes from the ENA db
#################################################################################

# install the needed packages
install.packages("RCurl", repos="http://cran.rstudio.com/")
install.packages("XML", repos="http://cran.rstudio.com/")
# SPARQL needs RCurl and XML
install.packages("SPARQL", repos="http://cran.rstudio.com/")

library("SPARQL")
library("XML")
library("RCurl")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
# genomes <- genome.and.organisms[,1]


###-------- Retrieving all gene sequences of genome ---------------###

# From the prodigal annotation tool we will retrieve all the annotated gene sequences
# We will go through all the genomes and all the gene sequences with its names 
# will be written to a file

query.gene.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?gene ?locus ?product 
WHERE {
VALUES ?sample { <http://gbol.life#xxx> }    
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:feature ?gene .
    ?gene gbol:transcript ?transcript .
    ?gene gbol:origin ?origin .
    ?origin <http://www.w3.org/ns/prov#provo:wasAttributedTo> ?att .
    ?att gbol:name "prodigal" .
    ?att gbol:version "2.6.3" .
    	?transcript gbol:feature ?feat .
    	?transcript gbol:locus_tag ?locus .
        ?transcript gbol:sequence ?CDS .
    	?transcript gbol:product ?product .
}    
'

ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

# creating the folder to save the data in
outfolder <- "CDS_data_withgenenames/"  
if (!file.exists(outfolder))dir.create(outfolder)

genomeID <- "GCA_000008865"

# takes every genome and writes the gene_ID and CDS to a file
for (genomeID in genome.and.organisms[,1]) { #always put he { on this line
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CDS_names.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) { 
    # substitute xxx with the genomeID in the domain.sparql query
    sub.gene.sparql <- sub("xxx", genomeID, query.gene.seqs)
    # running curl from command line
    curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.gene.sparql,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
    curl <- gsub(pattern = "\n", replacement = " ", x = curl)
    system(curl)
    output.genes <- rbind(sub.gene.sparql, read.csv("tmp.txt", sep = "\t"))
    if (length(output.genes[,1]) > 2){
    gene.data <- output.genes[-1,]
    gene.data.withoutcomma <- gsub(",", "" ,gene.data[,3])
    gene.data[,3] <- NULL
    gene.data$X.product <- gene.data.withoutcomma
    write.table(gene.data, file=fileout, append=F, sep = ",",
                row.names = F, quote=F, col.names=T )
    }
  }
}

