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

#---------Retrieving genomes from MicroDB--------#
genomes.Myc.Spir <- "
PREFIX ssb:<http://csb.wur.nl/genome/>
SELECT DISTINCT ?genome ?organism
WHERE {
  { ?genome a ssb:Genome .
  ?genome ssb:organism ?organism .           
	FILTER(regex(?organism, 'Spiroplasma', 'i'))
	}
UNION {
  ?genome a ssb:Genome .
  ?genome ssb:organism ?organism .
	FILTER(regex(?organism, 'Mycoplasma', 'i'))
  }
}
"

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"

# We want to retrieve the genomes with the location of it on the web
output <- SPARQL(url = endpoint, query = genomes.Myc.Spir)

# Edit genomenumbers, removing the hyperlink notation
genomes.and.organisms <- output$results
genome.number <- substr(genomes.and.organisms$genome, 27, 41)
genomes.and.organisms$genome <- genome.number


#write genomenumbers and organisms to file
write.table(genomes.and.organisms, file = "genomes_Mycoplasma&Spiroplasma.csv", 
            append = F, sep = ",", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)


