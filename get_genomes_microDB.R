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
genomes10 <- "
PREFIX ssb:<http://csb.wur.nl/genome/>
SELECT DISTINCT ?genome ?organism
WHERE {
?genome a ssb:Genome .
?genome ssb:organism ?organism .
} 
ORDER BY ?organism
LIMIT 10
"

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"

# We want to retrieve the genomes with the location of it on the web
output <- SPARQL(url = endpoint, query = genomes10)

# 
genomes10.with.names <- output$results

write.table(genomes10.with.names, file = "genomes10.csv", 
            append = F, sep = ",", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)


