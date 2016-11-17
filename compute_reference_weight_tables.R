#############################################################################
####  Author: Melanie van den Bosch
####  Script for retrieving ribosomal domains of a genome
####  and then computing a reference weight table for every genome
#################################################################################

# install the needed packages
install.packages("RCurl")
install.packages("XML")
# SPARQL needs RCurl and XML
install.packages("SPARQL")
# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")


# loading required library to execute SPARQL query
library("SPARQL")
# loading required library to compute the codon frequency
library("Biostrings")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes10.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomes10 <- genome.and.organisms[,1]


##-------- Computing reference weight table  ---------------###

# based on ribosomal protein domains
# We will go through all the genomes and the weight tables 
# will be written to a file
# creating a variable to store all the ribosomal proteins from all Acinetobacter baumannii genomes
# The output contains CDSsequence and proteinsequence, however 157 distinct CDS are obtained

ribosomal.seqs <- "
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
  #?protein ssb:sequence ?p_seq .
  ?cds ssb:tool 'prodigal' .
  ?protein ssb:feature ?domain .
  ?domain ssb:begin ?d_begin .
  ?domain ssb:end ?d_end .
  ?domain ssb:signature ?signature .
  ?signature biopax:xref ?xref .
  ?xref biopax:db 'pfam' .
  ?xref biopax:id ?Pfam_id .
  ?signature ssb:interpro_description ?domain_description
  FILTER(regex(?domain_description, 'ribosomal protein', 'i'))
}
"

endpoint <- "http://ssb2:9999/blazegraph/namespace/MicroDB/sparql/MicroDB/sparql"
#storing the output of the query.


outfolder <- "Reference_weight_tables/"  
if (!file.exists(outfolder))dir.create(outfolder)


for (genomeID in genomes10) { 
  fileout <- paste(outfolder, genomeID, ".csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    #substituting the genome numbers
    genome.sub <- sub("xxx", genomeID, ribosomal.seqs) 
    # Run SPARQL query for all genomes
    output.all <- SPARQL(url = endpoint, query = genome.sub) 
    #only retrieve results slice
    ribosomal.domain.data <- output.all$results
    # paste all the coding sequences together
    pasted.cdseqs <- paste(ribosomal.domain.data[,4], sep="", collapse="")
    # compute codon frequency
    codon.frequency <- trinucleotideFrequency(DNAString(pasted.cdseqs), as.prob=F, with.labels=T, as.array=F)
    
    
    #write.table(ribosomal.domain.data, file=fileout, append=F, sep = ",",
     #           row.names = F, quote=F, col.names=T ) # write the ribosomal domains of every genome to a file
  #remove the write.table, we want to convert the ribosomal.domain.data into a weight table
    }
}

output1 <- SPARQL(url = endpoint, query = AB1ribo, ns = c("ssb","http://csb.wur.nl/genome/"))
data1 <- data.frame(output1$results, stringsAsFactors=F)
# pastes all the sequences together
all_seq1 <-paste(as.matrix(data1)[,1], sep="", collapse="")
w_Abaum <- trinucleotideFrequency(DNAString(all_seq1), as.prob = F, with.labels = T, as.array = F)
str(w_Abaum)
w_Abaum

w_Abaumt <- t(w_Abaum)
# converting the weight vector to dataframe
w_data <- as.data.frame(w_Abaumt, stringsAsFactors = F)
