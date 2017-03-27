#!usr/bin/env Rscript

###################################################################
## Author: Melanie van den Bosch
## Script for computing the CAI of the catalase gene 
## and observing whether it is higher than the mean CAI of a genome
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

genomecount = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_cat.csv", sep="")
  #check if file already exists
  if (!file.exists(fileout)) {

    cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
      if (file.exists(cai.files)){
        cai.data <- read.csv(file = cai.files, header = TRUE, 
                           as.is=TRUE) #as.is to keep the it as char
        mean.cai <- mean(cai.data[,2])
      
        w.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
        if (file.exists(w.files)){
          w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
          # ordering on rownames, for it to be correctly used by cai function
          ordered.w <- w.data[with(w.data, order(V1)), ]
          # only leaving numbers
          w <- ordered.w[,2]
          

          
          
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
          catalase.seqs <- read.csv(file = fileout, sep = ",", header = TRUE, as.is = TRUE)
          # computing the postions of the domains with their corresponding postion in the CDS
          for (row in catalase.seqs){
            CDS_begin <- 3*(catalase.seqs[5]-1)+1
            CDS_end <- 3*(catalase.seqs[6])
          }
          
          # binding these columns to already existing dataframe
          # and converting to integers for later calculation
          catalase.seqs$CDS_domain_beginpos <- as.integer(CDS_begin[,1])
          catalase.seqs$CDS_domain_end <- as.integer(CDS_end[,1])
          
          # Substracting the CDS part that is associated to the protein domain
          # from columns 5 and 6 the positions of the begin and end of the 
          # CDS associated to protein domain are derived
          # when CDS are equal, take the integer from that column
          for (CDS in catalase.seqs[4]){
            dom_seq <- substr(CDS, catalase.seqs[catalase.seqs[,4] == CDS, 7], 
                              catalase.seqs[catalase.seqs[,4] == CDS, 8])
          }
          catalase.seqs$CDS_domain <- dom_seq 
          
          
          # shrinking the data to reduce computational time
          colnames(catalase.seqs) <- c("domain_ID", "domain_description", "genename", "CDS", 
                                       "d_begin","d_end", "CDS_beginpos", "CDS_end", "CDS_domain")
          keep <- c("domain_ID", "CDS_domain")
          catalase.seqs.CDS <- catalase.seqs[keep]
          
          cai.output <- compute.cai(catalase.seqs.CDS, genomeID, w, "tmpcatadomain.csv", "tmpcatadomain.fasta")
          catalase.cai <- mean(cai.output[,2])
          # draw the plot of mean cai and mean catalase cai
          xlim = c(0.1, 0.9)
          ylim = c(0.1, 0.9)
          col = "blue"
          if (genomecount == 0){  
            plot(catalase.cai, mean.cai, col=col, 
                 xlim = xlim, ylim = ylim,
                 pch = 18, main = "Average CAI vs. GC content",
                 xlab = "catalase CAI", 
                 ylab = "Mean CAI")
            grid(NULL, NULL, lty = 6, col = "cornsilk2")
            #legend("bottomright" ,c("dnaE1", "dnaE2/dnaE1", "polC/dnaE3", "Other"), cex=1.5, pch=18,
            #    col=c("red", "blue", "green", "grey") , bty="n")
            genomecount = genomecount + 1
          } else {
            points(catalase.cai, mean.cai, col=col, pch = 18)
          }
        }
      }
  }
}


