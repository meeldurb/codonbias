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
genome.and.organisms <- read.csv(file = "genomes10.csv", header = FALSE)
genomes10 <- genome.and.organisms[1]

#---------Retrieving domains from MicroDB--------#
# By going through all the genomes, we will retrieve the domains
# and all the data needed for further processing of the genomes
# these domains with all their information have to be written to a file




for (link in genomes10) 
{ 
  domain.sparql <-"
PREFIX ssb:<http://csb.wur.nl/genome/>
PREFIX biopax:<http://www.biopax.org/release/bp-level3.owl#>
SELECT DISTINCT ?Pfam_id ?d_begin ?d_end ?cds_seq ?p_seq
WHERE {
  VALUES ?genome { <hyperlink.genome> }
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
sub( "<hyperlink.genome>", toString(link), domain.sparql)
  
}

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"

# We want to retrieve the genomes with the location of it on the web
output <- SPARQL(url = endpoint, query = domain.sparql)

# 
domain.data <- output$results
