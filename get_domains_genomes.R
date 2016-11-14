##Author: Melanie van den Bosch
##Script for retrieving all the genomes from the 
##MicroDB namescpace on blazegraph


# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")

# loading required library to execute SPARQL query
library("SPARQL")

# reading the file with the genome names
genome.and.organisms <- read.csv(file = "genomes10.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomes10 <- genome.and.organisms[,1]

#---------Retrieving domains from MicroDB--------#
# By going through all the genomes, we will retrieve the domains
# and all the data needed for further processing of the genomes
# these domains with all their information have to be written to a file



domain.sparql <-"
PREFIX ssb:<http://csb.wur.nl/genome/>
PREFIX biopax:<http://www.biopax.org/release/bp-level3.owl#>
SELECT DISTINCT ?Pfam_id ?d_begin ?d_end ?cds_seq ?p_seq
WHERE {
  VALUES ?genome { <http://csb.wur.nl/genome/xxx> }
  ?genome a ssb:Genome .
  ?genome ssb:organism ?organism .
  ?genome ssb:dnaobject ?dna .
  ?dna ssb:feature ?cds .
  ?cds a ssb:Cds .
  ?cds ssb:sequence ?cds_seq .
  #?gene ssb:begin ?gene_begin .
  #?gene ssb:end ?gene_end .
  ?cds ssb:protein ?protein .
  ?protein ssb:sequence ?p_seq .
  ?cds ssb:tool 'prodigal' .
  ?protein ssb:feature ?Interpro .
  ?Interpro ssb:begin ?d_begin .
  ?Interpro ssb:end ?d_end .
  ?Interpro ssb:signature ?signatureAccession .
  ?signatureAccession biopax:xref ?xref .
  ?xref biopax:db 'pfam' .
  ?xref biopax:id ?Pfam_id .
}
LIMIT 10
"

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"


outfolder <- "Domain_data/"  
if (!file.exists(outfolder))dir.create(outfolder)

# takes every genome and writes the domain data to a file
for (genomeID in genomes10) { #always put he { on this line

 genome.sub <- sub("xxx", genomeID, domain.sparql)
  output.all <- SPARQL(url = endpoint, query = genome.sub)
 domain.data <- output.all$results
 #print (domain.data)
 #if (!file.exists(outfolder))dir.create(outfolder)
 fileout <- paste(outfolder, genomeID, sep="")
 #print (fileout)
 write.csv(domain.data, file=fileout)
 }


