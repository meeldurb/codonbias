########################################################
##Author: Melanie van den Bosch
##Script for retrieving all the domains of the genomes 
##from the ENA namespace on blazegraph
#########################################################

# install the needed packages
install.packages("RCurl")
install.packages("XML")
# load libraries
library("Rcurl")
library("XML")
library("SPARQL")

# set the working directory
setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomes <- genome.and.organisms[,1]

#---------Retrieving domains from MicroDB--------#
# By going through all the genomes, we will retrieve the domains
# and all the data needed for further processing of the genomes.
# The domains with all their information are written to a file



domain.sparql <-'
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?domain_id ?domain_description ?d_begin ?d_end ?CDS
WHERE {
VALUES ?sample { <http://gbol.life#GCA_000005825> }    
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:feature ?gene .
?gene gbol:transcript ?transcript .
?transcript gbol:sequence ?CDS .
?gene gbol:origin ?origin .
?origin <http://www.w3.org/ns/prov#provo:wasAttributedTo> ?annot .
?annot gbol:name "prodigal" .
?annot gbol:version "2.6.3".
?transcript gbol:feature ?cds .
?cds gbol:protein ?protein .
?protein gbol:feature ?feature .
?feature gbol:provenance ?provenance .
?provenance gbol:signature_description ?domain_description .
?feature gbol:location ?location .
?location gbol:begin ?beginiri .
?beginiri gbol:position ?d_begin .
?location gbol:end ?endiri .
?endiri gbol:position ?d_end .
?feature gbol:origin <http://gbol.life#/interproscan/interproscan-5.21-60.0> .
?feature gbol:xref ?xref .
?xref gbol:PrimaryAccession ?domain_id .
?xref gbol:db <http://identifiers.org/pfam/> .
} 
'

ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

# creating the folder to save the data in
outfolder <- "Domain_data_ENA_description/"  
if (!file.exists(outfolder))dir.create(outfolder)



# takes every genome and writes the Pfam ID, domain_begin, domain_end
# and CDS to a file
for (genomeID in genomes) { #always put he { on this line
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, ".csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) { 
    # substitute xxx with the genomeID in the domain.sparql query
    sub.domain.sparql <- sub("xxx", genomeID, domain.sparql)
    # running curl from command line
    curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.domain.sparql,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
    curl <- gsub(pattern = "\n", replacement = " ", x = curl)
    system(curl)
    output.domains <- rbind(sub.domain.sparql, read.csv("tmp.txt", sep = "\t"))
    domain.data <- output.domains[-1,]
    write.table(domain.data, file=fileout, append=F, sep = ",",
                row.names = F, quote=F, col.names=T )
  }
}


