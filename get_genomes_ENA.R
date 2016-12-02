##Author: Melanie van den Bosch
##Script for retrieving all the genomes from the 
##ENA database on blazegraph


# install the needed packages
install.packages("RCurl")
install.packages("XML")


#---------Retrieving genomes from MicroDB--------#
query.genomes <- "
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT  ?sample ?organism
WHERE {
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:scientificname ?organism .
}
LIMIT 10 
"

ENDPOINT = "http://ssb5.wurnet.nl:7200/repositories/ENA"

# We want to retrieve the genomes with the location of it on the web

curl = paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",query.genomes,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
curl = gsub(pattern = "\n", replacement = " ", x = curl)
system(curl)
output.genomes = rbind(query.genomes, read.csv("tmp.txt", sep = "\t"))

# omit the first row, which contains some NA values
output.genomes = output.genomes[-1,]
                

# Edit genomenumbers, removing the hyperlink notation
genome.number <- substr(output.genomes$X.sample, 19, 31)
output.genomes$X.sample <- genome.number


#write genomenumbers and organisms to file
write.table(output.genomes, file = "genomes_ENA10.csv", 
            append = F, sep = ",", row.names = FALSE, 
            quote = FALSE, col.names = FALSE)


