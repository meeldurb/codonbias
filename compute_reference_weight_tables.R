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
# loading required library to use CAI function
library("seqinr")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_orderedOngenome.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomes <- genome.and.organisms[,1]


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


amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                 "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                 "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                 "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                 "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Stop", "Cys", "Trp", "Cys",
                 "Leu", "Phe", "Leu", "Phe")

match.words = c("Mycoplasma", "Spiroplasma")

for (genomeID in genome.and.organisms[,1]) { #always put he { on this line
  # retrieves the index(i) of where the name is myco-/spiroplasma
  if (i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2] )) { 
    print (i)
  }
}

for (genomeID in genomes) { 
  fileout <- paste(outfolder, genomeID, ".csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    #substituting the genome numbers
    genome.sub <- sub("xxx", genomeID, ribosomal.seqs) 
    # Run SPARQL query for all genomes
    output.all <- SPARQL(url = endpoint, query = genome.sub)
    #only retrieve results slice
    ribosomal.domain.data <- output.all$results 
    #some do not have ribosomal domain data, should skip these
    if (length(ribosomal.domain.data) == 0) {
      next
    }
    # paste all the coding sequences together
    pasted.cdseqs <- paste(ribosomal.domain.data[,4], sep="", collapse="")
    # compute codon frequency
    codon.frequency <- trinucleotideFrequency(DNAString(pasted.cdseqs), as.prob=F, with.labels=T, as.array=F)
    # Overwriting all the codons
    # changing codons that code for the same amino acid and dividing by the sum of it
    # Phenylalanine
    codon.frequency[c("TTT", "TTC")]<- codon.frequency[c("TTT", "TTC")]/(sum(codon.frequency[c("TTT", "TTC")]))
    #Leucine
    codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]<- codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]/(sum(codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]))
    #Isoleucine
    codon.frequency[c("ATT", "ATC", "ATA")]<- codon.frequency[c("ATT", "ATC", "ATA")]/sum(codon.frequency[c("ATT", "ATC", "ATA")])
    #Valine
    codon.frequency[c("GTT", "GTC", "GTA", "GTG")]<- codon.frequency[c("GTT", "GTC", "GTA", "GTG")]/max(sum(codon.frequency[c("GTT", "GTC", "GTA", "GTG")]))
    #Serine
    codon.frequency[c("TCT", "TCC", "TCA", "TCG")]<- codon.frequency[c("TCT", "TCC", "TCA", "TCG")]/max(sum(codon.frequency[c("TCT", "TCC", "TCA", "TCG")]))
    #Proline
    codon.frequency[c("CCT", "CCC", "CCA", "CCG")]<- codon.frequency[c("CCT", "CCC", "CCA", "CCG")]/max(sum(codon.frequency[c("CCT", "CCC", "CCA", "CCG")]))
    #threonine
    codon.frequency[c("ACT", "ACC", "ACA", "ACG")]<- codon.frequency[c("ACT", "ACC", "ACA", "ACG")]/max(sum(codon.frequency[c("ACT", "ACC", "ACA", "ACG")]))
    #Alanine
    codon.frequency[c("GCT", "GCC", "GCA", "GCG")]<- codon.frequency[c("GCT", "GCC", "GCA", "GCG")]/max(sum(codon.frequency[c("GCT", "GCC", "GCA", "GCG")]))
    #Tyrosine
    codon.frequency[c("TAT", "TAC")]<- codon.frequency[c("TAT", "TAC")]/max(sum(codon.frequency[c("TAT", "TAC")]))
    #Histidine
    codon.frequency[c("CAT", "CAC")]<- codon.frequency[c("CAT", "CAC")]/max(sum(codon.frequency[c("CAT", "CAC")]))
    #glutamine
    codon.frequency[c("CAA", "CAG")]<- codon.frequency[c("CAA", "CAG")]/sum(codon.frequency[c("CAA", "CAG")])
    #asparginine
    codon.frequency[c("AAT", "AAC")]<- codon.frequency[c("AAT", "AAC")]/sum(codon.frequency[c("AAT", "AAC")])
    #lysine
    codon.frequency[c("AAA", "AAG")]<- codon.frequency[c("AAA", "AAG")]/sum(codon.frequency[c("AAA", "AAG")])
    #aspartic acid
    codon.frequency[c("GAT", "GAC")]<- codon.frequency[c("GAT", "GAC")]/sum(codon.frequency[c("GAT", "GAC")])
    #glutamic acid
    codon.frequency[c("GAA", "GAG")]<- codon.frequency[c("GAA", "GAG")]/sum(codon.frequency[c("GAA", "GAG")])
    #cysteine
    codon.frequency[c("TGT", "TGC")]<- codon.frequency[c("TGT", "TGC")]/sum(codon.frequency[c("TGT", "TGC")])
    #tryptophan
    codon.frequency["TGG"]<- codon.frequency["TGG"]/codon.frequency["TGG"]
    #arginine
    codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]<- codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]/sum(codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")])
    #serine
    codon.frequency[c("AGT", "AGC")]<- codon.frequency[c("AGT", "AGC")]/sum(codon.frequency[c("AGT", "AGC")])
    #glycine
    codon.frequency[c("GGT", "GGC", "GGA", "GGG")]<- codon.frequency[c("GGT", "GGC", "GGA", "GGG")]/sum(codon.frequency[c("GGT", "GGC", "GGA", "GGG")])
    #stopcodon
    codon.frequency[c("TAA", "TAG", "TGA")]<- codon.frequency[c("TAA", "TAG", "TGA")]/sum(codon.frequency[c("TAA", "TAG", "TGA")])
    #startcodon
    codon.frequency["ATG"]<- codon.frequency["ATG"]/codon.frequency["ATG"]
    
      #write.table(ribosomal.domain.data, file=fileout, append=F, sep = ",",
      #           row.names = F, quote=F, col.names=T ) # write the ribosomal domains of every genome to a file
      #remove the write.table, we want to convert the ribosomal.domain.data into a weight table
    codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
    codon.table <- cbind(codon.frequency.table, amino.acids)
    sorted.codon.table <- codon.table[order(amino.acids),]
    # correct the codon frequencies for the maximum codon.frequency belonging to one amino acid
    tableFinalC=NULL; 
    tableFinalV=NULL; 
    tableFinalA=NULL;
    for (aa in unique(sorted.codon.table[,2])) { 
      #cat(aa, "\n")
      which.max = which(sorted.codon.table[,2]==aa)
      tableFinalC = c(tableFinalC, rownames(sorted.codon.table)[which.max])
      tableFinalV = c(tableFinalV, sorted.codon.table[which.max,1]/max(sorted.codon.table[which.max,1] ))
      tableFinalA = c(tableFinalA, rep(aa, length(which.max)))
    }
    final.codon.table=data.frame(tableFinalC, tableFinalV, tableFinalA, stringsAsFactors = FALSE)
    colnames(final.codon.table) <- c("codon", "codon.frequency", "amino.acid")
    # write weight table to a file
    write.table(final.codon.table, file = fileout, append = F, sep = ",", row.names = F, quote = F, col.names = F)
  }
}



