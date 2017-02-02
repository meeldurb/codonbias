#############################################################################
####  Author: Melanie van den Bosch
####  Script for retrieving all genetic sequences of all genomes from the ENA db
#################################################################################

# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")

library("SPARQL")
library("XML")
library("RCurl")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
# genomes <- genome.and.organisms[,1]


###-------- Computing reference weight table  ---------------###

# based on ribosomal protein domains
# We will go through all the genomes and the weight tables 
# will be written to a file
# creating a variable to store all the ribosomal proteins from all genomes

query.gene.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?gene ?CDS 
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
        ?transcript gbol:sequence ?CDS .
}  
'

ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

# creating the folder to save the data in
outfolder <- "CDS_data/"  
if (!file.exists(outfolder))dir.create(outfolder)



# takes every genome and writes the Pfam ID, domain_begin, domain_end
# and CDS to a file
for (genomeID in genome.and.organisms[,1]) { #always put he { on this line
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_CDS.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) { 
    # substitute xxx with the genomeID in the domain.sparql query
    sub.gene.sparql <- sub("xxx", genomeID, query.gene.seqs)
    print(sub.gene.sparql)
    # running curl from command line
    curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.gene.sparql,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
    curl <- gsub(pattern = "\n", replacement = " ", x = curl)
    system(curl)
    output.genes <- rbind(sub.gene.sparql, read.csv("tmp.txt", sep = "\t"))
    gene.data <- output.genes[-1,]
    write.table(gene.data, file=fileout, append=F, sep = ",",
                row.names = F, quote=F, col.names=T )
  }
}

