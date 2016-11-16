##Author: Melanie van den Bosch
##Script for retrieving all the genomes from the 
##MicroDB namespace on blazegraph


# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")

# loading required library to execute SPARQL query
library("SPARQL")

# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_orderedOngenome.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomes <- genome.and.organisms[,1]

#---------Retrieving domains from MicroDB--------#
# By going through all the genomes, we will retrieve the domains
# and all the data needed for further processing of the genomes.
# The domains with all their information are written to a file



domain.sparql <-"
PREFIX ssb:<http://csb.wur.nl/genome/>
PREFIX biopax:<http://www.biopax.org/release/bp-level3.owl#>
SELECT DISTINCT ?Pfam_id ?d_begin ?d_end ?cds_seq
WHERE {
VALUES ?genome { <http://csb.wur.nl/genome/xxx> }
?genome a ssb:Genome .
?genome ssb:organism ?organism .
?genome ssb:dnaobject ?dna .
?dna ssb:feature ?cds .
?cds a ssb:Cds .
?cds ssb:sequence ?cds_seq .
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
"

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"

# creating the folder to save the data in
outfolder <- "Domain_data/"  
if (!file.exists(outfolder))dir.create(outfolder)

# takes every genome and writes the Pfam ID, domain_begin, domain_end
# and CDS to a file
for (genomeID in genomes) { #always put he { on this line
  fileout <- paste(outfolder, genomeID, ".csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    genome.sub <- sub("xxx", genomeID, domain.sparql)
    output.all <- SPARQL(url = endpoint, query = genome.sub)
    domain.data <- output.all$results
    write.table(domain.data, file=fileout, append=F, sep = ",",
              row.names = F, quote=F, col.names=T )
  }
}


