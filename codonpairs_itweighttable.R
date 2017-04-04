#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables of codon pairs 
####  Purpose of script: Computation of iterative weight tables 
####  of codon pairs
#################################################################################

# install needed packages
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")
# for strings
install.packages("stringr", repos="http://cran.rstudio.com/")

# loading required libraries
library("Biostrings")
library("seqinr")
library("stringr")

#loading my own functions
#source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")
#source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/ccodpairsweight.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAIpairs.R")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")




# open and read file with genomeIDs
genomeID <- "GCA_000003925"
genome.and.organisms <- read.csv(file = "test_genomes_ENA10.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


outfolderw <- "codonpairs_weight/"  
if (!file.exists(outfolderw))dir.create(outfolderw)

outfolderitw <- "codonpairs_itweight/"
if (!file.exists(outfolderitw))dir.create(outfolderitw)

outfolder100 <- "top100genesit_results/"
if (!file.exists(outfolder100))dir.create(outfolder100)

# looping through all the genomes
for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  w.it.files <- paste(outfolderitw, genomeID, "_witcodpairs.csv", sep="")
  w.files <- paste(outfolderw, genomeID, "_wcodpairs.csv", sep="")
  # if the initial w table is not written yet, do so
  if (!file.exists(w.files)) {
    #####~~~~~~~~~~~~~~~~~~~~~~~~~ Retrieving ribosomal protein genes ~~~~~~~~~~~~~~~~~~~~~~~~~#####
    
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
      print ("no ribosomal seqs found")
    }
    print ("new w table will be written")
    w.codpairtable <- compute.codpairs.weight(ribosomal.seqs.data[,4], genomeID)
    write.table(w.codpairtable, file = w.files, append = F, sep = ",", row.names = F, quote = F, col.names = F)
  } 
  # if the initial w table exists, open it
  if (!file.exists(w.it.files)) {
    
  if(file.exists(w.files)) {
    gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
    if (file.exists(gene.files)){
      w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
      gene.data <- read.csv(file = gene.files, header = TRUE, 
                            as.is=TRUE) #as.is to keep the it as char 
      cai.ini <- compute.codpairs.cai(gene.data, w.data, genomeID)
      
      # sort on CAI value and take top 100
      #sort.cai <- cai.data[order(-cai.data[,2]),]
      ini.sort.cai <- cai.ini[order(-cai.ini[,2]),]
      ini.top100 <- head(ini.sort.cai, 100)
      
      # We retrieve only the CDS of the gene_IDs that are in the top 100
      match.id <- as.vector(ini.top100[,1])
      gene.match <- gene.data[gene.data[,1] %in% match.id, ] 
      
      # then re-compute weight tables 
      new.w.table <- compute.codpairs.weight(gene.match[,2], genomeID)

      # We re-calculate the CAI of all the gene_IDs 
      cai.res <- compute.codpairs.cai(gene.data, new.w.table, genomeID)
      
      # and take the top 100 again
      res.sort.cai <- cai.res[order(-cai.res[,2]),] 
      res.top100 <- head(res.sort.cai, 100)
      
      #compare initial and result top 100
      diff.count <- length(setdiff(ini.top100[,1], res.top100[,1]))
      cat(paste("differences between tables is ", diff.count, "\n"))
      
      # if the difference between both lists is lower than 5, run the analysis again
      # else save the resulting weight table
      it.count = 1
      while (diff.count > 0 && it.count <= 20){
        ini.top100 <- res.top100
        # We retrieve only the CDS of the gene_IDs that are in the top 100
        match.id <- as.vector(ini.top100[,1])
        gene.match <- gene.data[gene.data[,1] %in% match.id, ] 
        
        # then re-compute weight tables 
        new.w.table <- compute.codpairs.weight(gene.match[,2], genomeID)

        
        # We re-calculate the CAI of all the gene_IDs 
        cai.res <- compute.codpairs.cai(gene.data, new.w.table, genomeID)
        
        # and take the top 100 again
        res.sort.cai <- cai.res[order(-cai.res[,2]),] 
        res.top100 <- head(res.sort.cai, 100)
        
        #compare initial and result top 100
        diff.count <- length(setdiff(ini.top100[,1], res.top100[,1]))
        cat(paste("differences between tables is ", diff.count, "\n"))
        
        # keep a count of the iterations, loop needs to stop after 20
        it.count <- sum(it.count, 1)
      }

      # write weight table to file
      write.table(new.w.table, file = w.it.files, append = FALSE, sep = ",", 
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
      write.table(res.top100, file = paste(outfolder100, genomeID, "_restop100.csv", sep=""),
                  append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      }
    }
  }
}





#####___________________ example script ____________________####

#############______________________ fake sequences ______________________###########
# # counting the codon pairs
# # fake sequence
# seq1 <- "GGTGGTGGTGGCGGTGGCGGTGGCGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGG"
# seq2 <- "GGAGGAGGAGGAGGGGGAGGGGGAGGGGGAGGAGGTGGAGGTGGAGGTGGAGGTGGAGGT"
# seq3 <- "GGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGG"
# seq4 <- "AACTGCAACTGCAACTGCAACTGCTACAGGTACAGGTACAGGTACAGGTACAGGTACAGG"
# # how it should be
# #seqs3 <- c(seq1, seq2)
# # convert to df, because also read like this when loading the data
# seqs <- as.data.frame(c(seq1, seq2, seq3, seq4))
# seqs2 <- as.vector(seqs[,1])
#############______________________ fake sequences ______________________###########


# getting the files
genomeID <- "GCA_000003925"
gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
gene.data <- read.csv(file = gene.files, header = TRUE, 
                      as.is=TRUE) #as.is to keep the it as char


# adding extra column to df for the counting of codonpairs
counter <- rep(0, length(codonpair.table[,1]))
codonpair.table$count <- counter

seqcount <- 1
for (DNAseq in gene.data[,2]){
  # after every codon place an "_"
  sq <- gsub("(.{3})", "\\1 ", DNAseq)
  sq <- sub("\\s+$", "", sq)
  sq <- gsub(" ", "_", sq)
  seqcount <- sum(seqcount, 1)
  print (seqcount)
  for (codonpair in codonpair.table[,2]){
    pair.count <- str_count(sq, pattern=codonpair)
    if (pair.count > 0){
      # take the codonpair from data that is the same as this one
      sel <- which(codonpair == codonpair.table[,2])
      codonpair.table[sel,3] <- sum(codonpair.table[sel,3], pair.count)

    }
  }
}
codonpair.table$frequency <- rep(0, length(codonpair.table[,1]))
codonpair.table <- codonpair.table[order(codonpair.table[,1]), ]
unique.aa <- unique(codonpair.table[,1])
codonrow <- 1
for (aa in unique.aa){
  aapair.count <- codonpair.table[which(codonpair.table[,1] == aa), ]
  totcount.aa <- sum(aapair.count[,3])
  print (aa)
  #print (totcount.aa)
  for (count in aapair.count[,3]) {
    print (codonpair.table[meh,2])
    print (count)
    freq <- count/totcount.aa
    print (freq)
    codonpair.table[codonrow, 4] <- freq
    meh <- sum(codonrow, 1)
  }
}






# check whether codons are saved in dataframe
findcod <- which("AAC_TGC" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("TGC_AAC" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("TAC_AGG" == codonpair.table[,2])
codonpair.table[findcod,3]
findcod <- which("AGG_TAC" == codonpair.table[,2])
codonpair.table[findcod,3]


#fakedata
data <- sample(1:100, 25)
names(data) <- pairsCodons
data
AAp <- unique(tranlPairs)

for (aap in AAp){
  sel <- names(which(tranlPairs==aap))
  data[sel]
}

AAp
pairsCodons
data

