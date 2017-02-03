#############################################################################
####  Author: Melanie van den Bosch
####  Script for retrieving ribosomal domains of a genome for ENA db
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

# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
# genomes <- genome.and.organisms[,1]


###-------- Computing reference weight table  ---------------###

# based on ribosomal protein domains
# We will go through all the genomes and the weight tables 
# will be written to a file
# creating a variable to store all the ribosomal proteins from all genomes

query.ribosomal.seqs <- '
PREFIX gbol: <http://gbol.life#>
SELECT DISTINCT ?domain_id ?d_begin ?d_end ?CDS
WHERE {
VALUES ?sample { <http://gbol.life#xxx> }    
?sample a gbol:Sample .
?contig gbol:sample ?sample .
?contig gbol:feature ?gene .
?gene gbol:transcript ?transcript .
?transcript gbol:sequence ?CDS .
?gene gbol:origin ?origin .
?origin <http://www.w3.org/ns/prov#provo:wasAttributedTo> ?annotation .
?annotation gbol:name "prodigal" .
?annotation gbol:version "2.6.3".
?transcript gbol:feature ?cds .
?cds gbol:protein ?protein .
?protein gbol:feature ?feature .
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
FILTER(regex(?domain_description, "ribosomal protein", "i"))
}
'

ENDPOINT = "http://ssb2.wurnet.nl:7201/repositories/ENA"

# creating the folder to save the data in
outfolder <- "Reference_weight_tables_ENA/"  
if (!file.exists(outfolder))dir.create(outfolder)


# amino acid names of standard genetic code and myc/spiroplasma genetic code
amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                 "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                 "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                 "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                 "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Stop", "Cys", "Trp", "Cys",
                 "Leu", "Phe", "Leu", "Phe")

myc.spir.amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                 "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                 "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                 "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                 "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Trp", "Cys", "Trp", "Cys",
                 "Leu", "Phe", "Leu", "Phe")


# # takes every genome and retrieves the ribosomal sequences.
## then computes the codon frequency of all the codons in these sequences
## writes the weight table of all genomes seperately to a file

for (genomeID in genome.and.organisms[,1]) { 
  fileout <- paste(outfolder, genomeID, ".csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {
    # substitute xxx with the genomeID in the query.ribosomal.seqs query
    sub.query.ribosomal.seqs <- sub("xxx", genomeID, query.ribosomal.seqs)
    # running curl from command line
    curl <- paste0("curl -s -X POST ",ENDPOINT," --data-urlencode 'query=",sub.query.ribosomal.seqs,"' -H 'Accept:text/tab-separated-values' > tmp.txt")
    curl <- gsub(pattern = "\n", replacement = " ", x = curl)
    system(curl)
    output.riboseqs <- rbind(sub.query.ribosomal.seqs, read.csv("tmp.txt", sep = "\t"))
    # slice off 1st row which contains NA values
    ribosomal.seqs.data <- output.riboseqs[-1,]
    #some genomes do not have ribosomal domain data, we should skip these
    if (length(ribosomal.seqs.data) == 0) {
      next
    }
    # paste all the coding sequences together
    pasted.cdseqs <- paste(ribosomal.seqs.data[,4], sep="", collapse="")
    # compute codon frequency
    codon.frequency <- trinucleotideFrequency(DNAString(pasted.cdseqs), as.prob=F, with.labels=T, as.array=F)
   
    # Overwriting all the codons
    # changing codons that code for the same amino acid and dividing by the sum of it
    # Phenylalanine
    codon.frequency[c("TTT", "TTC")]<- codon.frequency[c("TTT", "TTC")]/sum(codon.frequency[c("TTT", "TTC")])
    #Leucine
    codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]<- codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]/sum(codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")])
    #Isoleucine
    codon.frequency[c("ATT", "ATC", "ATA")]<- codon.frequency[c("ATT", "ATC", "ATA")]/sum(codon.frequency[c("ATT", "ATC", "ATA")])
    #Valine
    codon.frequency[c("GTT", "GTC", "GTA", "GTG")]<- codon.frequency[c("GTT", "GTC", "GTA", "GTG")]/sum(codon.frequency[c("GTT", "GTC", "GTA", "GTG")])
    #Serine
    codon.frequency[c("TCT", "TCC", "TCA", "TCG")]<- codon.frequency[c("TCT", "TCC", "TCA", "TCG")]/sum(codon.frequency[c("TCT", "TCC", "TCA", "TCG")])
    #Proline
    codon.frequency[c("CCT", "CCC", "CCA", "CCG")]<- codon.frequency[c("CCT", "CCC", "CCA", "CCG")]/sum(codon.frequency[c("CCT", "CCC", "CCA", "CCG")])
    #threonine
    codon.frequency[c("ACT", "ACC", "ACA", "ACG")]<- codon.frequency[c("ACT", "ACC", "ACA", "ACG")]/sum(codon.frequency[c("ACT", "ACC", "ACA", "ACG")])
    #Alanine
    codon.frequency[c("GCT", "GCC", "GCA", "GCG")]<- codon.frequency[c("GCT", "GCC", "GCA", "GCG")]/sum(codon.frequency[c("GCT", "GCC", "GCA", "GCG")])
    #Tyrosine
    codon.frequency[c("TAT", "TAC")]<- codon.frequency[c("TAT", "TAC")]/sum(codon.frequency[c("TAT", "TAC")])
    #Histidine
    codon.frequency[c("CAT", "CAC")]<- codon.frequency[c("CAT", "CAC")]/sum(codon.frequency[c("CAT", "CAC")])
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
    #arginine
    codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]<- codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]/sum(codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")])
    #serine
    codon.frequency[c("AGT", "AGC")]<- codon.frequency[c("AGT", "AGC")]/sum(codon.frequency[c("AGT", "AGC")])
    #glycine
    codon.frequency[c("GGT", "GGC", "GGA", "GGG")]<- codon.frequency[c("GGT", "GGC", "GGA", "GGG")]/sum(codon.frequency[c("GGT", "GGC", "GGA", "GGG")])
    #startcodon
    codon.frequency["ATG"]<- codon.frequency["ATG"]/codon.frequency["ATG"]
    
    # Searching for Mycoplasma and Spiroplasma, using other weight tables
    match.words <- c("Mycoplasma", "Spiroplasma")
    # i contains the indices where Myco/Spiro is found
    i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
    genomesMycoSpiro <- genome.and.organisms[i,1]
    # when the genome number of myc/spir is not found it will run the default
    # else it will run the codon table of myc+spiroplasma
    if(!(genomeID %in% genomesMycoSpiro)) {
       #tryptophan (default)
        codon.frequency["TGG"]<- codon.frequency["TGG"]/codon.frequency["TGG"]
        #stopcodon (default)
        codon.frequency[c("TAA", "TAG", "TGA")]<- codon.frequency[c("TAA", "TAG", "TGA")]/sum(codon.frequency[c("TAA", "TAG", "TGA")])
        codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
        codon.table <- cbind(codon.frequency.table, amino.acids)
        sorted.codon.table <- codon.table[order(amino.acids),]
        } else { 
        #tryptophan (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
        codon.frequency[c("TGG", "TGA")]<- codon.frequency[c("TGG", "TGA")]/sum(codon.frequency[c("TGG", "TGA")])
        #stopcodon (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
        codon.frequency[c("TAA", "TAG")]<- codon.frequency[c("TAA", "TAG")]/sum(codon.frequency[c("TAA", "TAG")])
        codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
        codon.table <- cbind(codon.frequency.table, myc.spir.amino.acids)
        sorted.codon.table <- codon.table[order(myc.spir.amino.acids),]
        }
    
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


#############_________________Function for computing refweighttables_________________#############

compute.weight <- function(seqs.df, genome_ID){

  pasted.cdseqs <- paste(seqs.df, sep="", collapse="")
  # compute codon frequency
  codon.frequency <- trinucleotideFrequency(DNAString(pasted.cdseqs), as.prob=F, with.labels=T, as.array=F)
  
  
  # amino acid names of standard genetic code and myc/spiroplasma genetic code
  amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                   "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                   "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                   "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                   "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Stop", "Cys", "Trp", "Cys",
                   "Leu", "Phe", "Leu", "Phe")
  
  myc.spir.amino.acids <- c("Lys", "Asn", "Lys", "Asn", "Thr", "Thr", "Thr", "Thr", "Arg", "Ser", "Arg", "Ser", 
                            "Ile", "Ile", "Met/Start", "Ile", "Gln", "His", "Gln", "His", "Pro", "Pro", "Pro", "Pro",
                            "Arg", "Arg", "Arg", "Arg", "Leu", "Leu", "Leu", "Leu", "Glu", "Asp", "Glu", "Asp",
                            "Ala", "Ala", "Ala", "Ala", "Gly", "Gly", "Gly", "Gly", "Val", "Val", "Val", "Val",
                            "Stop", "Tyr", "Stop", "Tyr", "Ser", "Ser", "Ser", "Ser", "Trp", "Cys", "Trp", "Cys",
                            "Leu", "Phe", "Leu", "Phe")
  
  
  # Overwriting all the codons
  # changing codons that code for the same amino acid and dividing by the sum of it
  # Phenylalanine
  codon.frequency[c("TTT", "TTC")]<- codon.frequency[c("TTT", "TTC")]/sum(codon.frequency[c("TTT", "TTC")])
  #Leucine
  codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]<- codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")]/sum(codon.frequency[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")])
  #Isoleucine
  codon.frequency[c("ATT", "ATC", "ATA")]<- codon.frequency[c("ATT", "ATC", "ATA")]/sum(codon.frequency[c("ATT", "ATC", "ATA")])
  #Valine
  codon.frequency[c("GTT", "GTC", "GTA", "GTG")]<- codon.frequency[c("GTT", "GTC", "GTA", "GTG")]/sum(codon.frequency[c("GTT", "GTC", "GTA", "GTG")])
  #Serine
  codon.frequency[c("TCT", "TCC", "TCA", "TCG")]<- codon.frequency[c("TCT", "TCC", "TCA", "TCG")]/sum(codon.frequency[c("TCT", "TCC", "TCA", "TCG")])
  #Proline
  codon.frequency[c("CCT", "CCC", "CCA", "CCG")]<- codon.frequency[c("CCT", "CCC", "CCA", "CCG")]/sum(codon.frequency[c("CCT", "CCC", "CCA", "CCG")])
  #threonine
  codon.frequency[c("ACT", "ACC", "ACA", "ACG")]<- codon.frequency[c("ACT", "ACC", "ACA", "ACG")]/sum(codon.frequency[c("ACT", "ACC", "ACA", "ACG")])
  #Alanine
  codon.frequency[c("GCT", "GCC", "GCA", "GCG")]<- codon.frequency[c("GCT", "GCC", "GCA", "GCG")]/sum(codon.frequency[c("GCT", "GCC", "GCA", "GCG")])
  #Tyrosine
  codon.frequency[c("TAT", "TAC")]<- codon.frequency[c("TAT", "TAC")]/sum(codon.frequency[c("TAT", "TAC")])
  #Histidine
  codon.frequency[c("CAT", "CAC")]<- codon.frequency[c("CAT", "CAC")]/sum(codon.frequency[c("CAT", "CAC")])
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
  #arginine
  codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]<- codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")]/sum(codon.frequency[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")])
  #serine
  codon.frequency[c("AGT", "AGC")]<- codon.frequency[c("AGT", "AGC")]/sum(codon.frequency[c("AGT", "AGC")])
  #glycine
  codon.frequency[c("GGT", "GGC", "GGA", "GGG")]<- codon.frequency[c("GGT", "GGC", "GGA", "GGG")]/sum(codon.frequency[c("GGT", "GGC", "GGA", "GGG")])
  #startcodon
  codon.frequency["ATG"]<- codon.frequency["ATG"]/codon.frequency["ATG"]
  
  #load all the genomeIDs and organisms for checking for Myco/Spiroplasma
  genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                   as.is=TRUE) #as.is to keep the it as char
  
  # Searching for Mycoplasma and Spiroplasma, which use other genetic codon table
  match.words <- c("Mycoplasma", "Spiroplasma")
  # i contains the indices where Myco/Spiro is found
  i <- grep(paste(match.words, collapse="|"), genome.and.organisms[,2])
  genomesMycoSpiro <- genome.and.organisms[i,1]
  # when the genome number of myc/spir is not found it will run the standard genetic code
  # else it will run the codon table of myc+spiroplasma
  if(!(genomeID %in% genomesMycoSpiro)) {
    #tryptophan (default)
    codon.frequency["TGG"]<- codon.frequency["TGG"]/codon.frequency["TGG"]
    #stopcodon (default)
    codon.frequency[c("TAA", "TAG", "TGA")]<- codon.frequency[c("TAA", "TAG", "TGA")]/sum(codon.frequency[c("TAA", "TAG", "TGA")])
    codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
    codon.table <- cbind(codon.frequency.table, amino.acids)
    sorted.codon.table <- codon.table[order(amino.acids),]
  } else { 
    #tryptophan (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
    codon.frequency[c("TGG", "TGA")]<- codon.frequency[c("TGG", "TGA")]/sum(codon.frequency[c("TGG", "TGA")])
    #stopcodon (Myco+Spiro; the stop codon TGA is a W in mycoplasma)
    codon.frequency[c("TAA", "TAG")]<- codon.frequency[c("TAA", "TAG")]/sum(codon.frequency[c("TAA", "TAG")])
    codon.frequency.table <- as.data.frame(codon.frequency, stringsAsFactors = F)
    codon.table <- cbind(codon.frequency.table, myc.spir.amino.acids)
    sorted.codon.table <- codon.table[order(myc.spir.amino.acids),]
  }
  
  # correct the codon frequencies for the maximum codon.frequency belonging to one amino acid
  tableFinalC=NULL; 
  tableFinalV=NULL; 
  tableFinalA=NULL;
  for (aa in unique(sorted.codon.table[,2])) { 
    #cat(aa, "\n")
    which.max <- which(sorted.codon.table[,2]==aa)
    tableFinalC <- c(tableFinalC, rownames(sorted.codon.table)[which.max])
    tableFinalV <- c(tableFinalV, sorted.codon.table[which.max,1]/max(sorted.codon.table[which.max,1] ))
    tableFinalA <- c(tableFinalA, rep(aa, length(which.max)))
  }
  final.codon.table <- data.frame(tableFinalC, tableFinalV, tableFinalA, stringsAsFactors = FALSE)
  colnames(final.codon.table) <- c("codon", "codon.frequency", "amino.acid")
  return(final.codon.table)
}

codon.table <- compute.weight(ribosomal.seqs.data[,4], "GCA_000003645")
