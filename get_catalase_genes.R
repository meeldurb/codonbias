#!usr/bin/env Rscript

###################################################################
## Author: Melanie van den Bosch
## Script for retrieval of the catalase genes from the ENA db
###################################################################

# Packages that need to be installed
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")

# loading required libraries
library("Biostrings")
library("seqinr")

source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")



genomeID <- "GCA_000005825"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

outfolder <- "Catalase_genes/"  
if (!file.exists(outfolder))dir.create(outfolder)


#####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving catalase protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####

query.catalase.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?domain_id ?domain_description ?genename ?CDS ?d_begin ?d_end
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
  fileout <- paste(outfolder, genomeID, "_cat.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    # check whether CAI and w files exists, else the catalase CAI cannot be calculated anyway
    cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
      if (file.exists(cai.files)){
       w.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
        if (file.exists(w.files)){
          #####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving catalase protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####
          
          sub.query.catalase.seqs <- sub("xxx", genomeID, query.catalase.seqs)
          # running curl from command line
          curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.query.catalase.seqs,
                         "' -H 'Accept:text/tab-separated-values' > tmp.txt")
          curl <- gsub(pattern = "\n", replacement = " ", x = curl)
          system(curl)
          output.catseqs <- rbind(sub.query.catalase.seqs, read.csv("tmp.txt", sep = "\t"))
          # slice off 1st row which contains NA values
          catalase.seqs.data <- output.catseqs[-1,]
          #some genomes do not have catalase domain data, we should skip these
          if (length(catalase.seqs.data) == 0) {
            cat (genomeID, "has no catalase genes\n")
            next 
          }
          write.table(catalase.seqs.data, file = fileout, append = F, sep = ",", row.names = F, quote = F, col.names = T)
        }
      }
  }
}
  
  

  
  
  
