#!usr/bin/env Rscript

###################################################################
## Author: Melanie van den Bosch
## Script for computing the CAI of the catalase gene 
## and observing whether it is higher than the mean CAI of a genome
###################################################################



genomeID <- "GCA_000005825"


#####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving catalase protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####

query.catalase.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?genename ?CDS ?domain_id ?domain_description ?d_begin ?d_end
WHERE {
VALUES ?sample { <http://gbol.life#xxx> }    
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:feature ?gene .
?gene gbol:transcript ?transcript .
?transcript gbol:sequence ?CDS .
?transcript gbol:gene ?genename  .
?gene gbol:origin ?origin .
?origin <http://www.w3.org/ns/prov#provo:wasAttributedTo> ?annotation .
?annotation gbol:name "prodigal" .
?annotation gbol:version "2.6.3".
?transcript gbol:feature ?cds .
?cds gbol:protein ?protein .
?protein gbol:sequence ?pseq .
?protein gbol:feature ?feature .
?protein gbol:xref ?xref .
?xref gbol:PrimaryAccession ?accession.
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
FILTER(regex(?domain_description, "Catalase"))
FILTER(!regex(?domain_description, " "))
}
'
  
ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
    if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, header = TRUE, 
                          as.is=TRUE) #as.is to keep the it as char
    mean.cai <- mean(cai.data[,2])
   
    #####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving catalase protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####
    
      sub.query.catalase.seqs <- sub("xxx", genomeID, query.catalase.seqs)
      # running curl from command line
      curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.query.catalase.seqs,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
      curl <- gsub(pattern = "\n", replacement = " ", x = curl)
      system(curl)
      output.catseqs <- rbind(sub.query.catalase.seqs, read.csv("tmp.txt", sep = "\t"))
      # slice off 1st row which contains NA values
      catalase.seqs.data <- output.catseqs[-1,]
      #some genomes do not have catalase domain data, we should skip these
      if (length(catalase.seqs.data) == 0) {
        next 
      }
    }
  } 

