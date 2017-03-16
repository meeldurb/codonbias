#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables 
####  Purpose of script: Re-computation of weight tables by selecting the top 25
####  domains with highest cai after each iteration is checked whether the top 25 
####  domains are comparable to previous iteration
#################################################################################

# install needed packages
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")

# loading required libraries
library("Biostrings")
library("seqinr")

#loading my own functions
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cweight.R")
source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

#genomeID <- "GCA_000003925"

# creating the folder to save the data in
outfolder <- "robustnesscheck_iterated_weight_tables/"  
if (!file.exists(outfolder))dir.create(outfolder)

genomeID.table <- NULL
itcount.table <- NULL
diffcount.table <- NULL

for (genomeID in genome.and.organisms[,1]) { 
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_robust_it.csv", sep="")
  if (!file.exists(fileout)) {

      #CDS data
      #cai.files <- paste("CAI_CDS/", genomeID, "_CAI_CDS.csv", sep = "")
      #cai.data <- read.csv(file = cai.files, header = TRUE, as.is = TRUE)
      #gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
      gene.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
      gene.data <- read.csv(file = gene.files, header = TRUE, 
                            as.is=TRUE) #as.is to keep the it as char
      
      #compute initial weight table from 25 randomly selected genes
      randomgenes <- sample(gene.data[,2], 25)
      w.data <- compute.weight(randomgenes, genomeID)
          # order on codon because of cai function
          ordered.w <- w.data[with(w.data, order(w.data[,1])), ]
          # only leaving numbers
          w <- ordered.w[,2]
          
      # compute CAI for the first round
      cai.ini <- compute.cai(gene.data, genomeID, w, "tmpc.csv", "tmpc.fasta")

      # sort on CAI value and take top 25
      #sort.cai <- cai.data[order(-cai.data[,2]),]
      ini.sort.cai <- cai.ini[order(-cai.ini[,2]),]
      ini.top25 <- head(ini.sort.cai, 25)
      
      # We retrieve only the CDS of the gene_IDs that are in the top 25
      match.id <- as.vector(ini.top25[,1])
      gene.match <- gene.data[gene.data[,1] %in% match.id, ] 
      
      # then re-compute weight tables 
      w.table <- compute.weight(gene.match[,2], genomeID)
      ordered.w <- w.table[with(w.table, order(w.table[,1])), ]
      # only leaving numbers
      w <- ordered.w[,2]
      
      
      # We re-calculate the CAI of all the gene_IDs 
      cai.res <- compute.cai(gene.data, genomeID, w, "tmpc.csv", "tmpc.fasta")
      
      # and take the top 25 again
      res.sort.cai <- cai.res[order(-cai.res[,2]),] 
      res.top25 <- head(res.sort.cai, 25)
      
      #compare initial and result top 25
      diff.count <- length(setdiff(ini.top25[,1], res.top25[,1]))
      cat(paste("differences between tables is ", diff.count, "\n"))
      
      # if the difference between both lists is lower than 5, run the analysis again
      # else save the resulting weight table
      it.count = 1
      while (diff.count > 0 && it.count <= 20){
        ini.top25 <- res.top25
        # We retrieve only the CDS of the gene_IDs that are in the top 25
        match.id <- as.vector(ini.top25[,1])
        gene.match <- gene.data[gene.data[,1] %in% match.id, ] 
        
        # then re-compute weight tables 
        w.table <- compute.weight(gene.match[,2], genomeID)
        ordered.w <- w.table[with(w.table, order(w.table[,1])), ]
        # only leaving numbers
        w <- ordered.w[,2]
        
        # We re-calculate the CAI of all the gene_IDs 
        cai.res <- compute.cai(gene.data, genomeID, w, "tmpc.csv", "tmpc.fasta")
        
        # and take the top 25 again
        res.sort.cai <- cai.res[order(-cai.res[,2]),] 
        res.top25 <- head(res.sort.cai, 25)
        
        #compare initial and result top 25
        diff.count <- length(setdiff(ini.top25[,1], res.top25[,1]))
        cat(paste("differences between tables is ", diff.count, "\n"))
        
        # keep a count of the iterations, loop needs to stop after 20
        it.count <- sum(it.count, 1)
      }
      # save count of iterations for each genome
      genomeID.table <- c(genomeID.table, genomeID)
      itcount.table <- c(itcount.table, it.count)
      diffcount.table <- c(diffcount.table, diff.count)
      
      
      # write weight table to file
      write.table(w.table, file = fileout, append = FALSE, sep = ",", 
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
      write.table(res.top25, file = paste(outfolder25, genomeID, "_restop25.csv", sep=""),
                  append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      }
    }
  }
}

# fill dataframe of iteration count per genome after all weight tables were computed
iteration.df <- data.frame(genomeID = genomeID.table, iterations = itcount.table, 
                           difference = diffcount.table, stringsAsFactors = FALSE)

write.table(iteration.df, file = "iterationcount.csv", append = FALSE, sep = ",", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)






