##Author: Melanie van den Bosch
##Script for retrieving all the genomes from the 
##MicroDB namescpace on blazegraph


# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")



#---------Retrieving genomes from MicroDB--------#
query.genomes <- "
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT  ?sample ?organism
WHERE {
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:scientificname ?organism .
}
"

ENDPOINT = "http://ssb5.wurnet.nl:7200/repositories/ENA"

# We want to retrieve the genomes with the location of it on the web

curl = paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",query.genomes,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
curl = gsub(pattern = "\n", replacement = " ", x = curl)
system(curl)
output.genomes = rbind(query.genomes, read.csv("tmp.txt", sep = "\t"))

#omit the first row
output.genomes = output.genomes[-1,]
                



# Edit genomenumbers, removing the hyperlink notation
genomes.and.organisms <- output$results
genome.number <- substr(genomes.and.organisms$genome, 27, 41)
genomes.and.organisms$genome <- genome.number


#write genomenumbers and organisms to file
write.table(genomes.and.organisms, file = "genomes_orderedOngenomes.csv", 
            append = F, sep = ",", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)


