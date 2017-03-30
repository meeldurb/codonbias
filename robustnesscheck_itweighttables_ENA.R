#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Computing iterated reference weight tables and checking for robustness 
####  Purpose of script: Re-computation of weight tables by selecting 25 random
####  selected genes as seed. Then computing the weight tables by iterative
####  algorithm
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
      #write.table(res.top25, file = paste(outfolder25, genomeID, "_restop25.csv", sep=""),
       #           append = FALSE, sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
      }
}


# # fill dataframe of iteration count per genome after all weight tables were computed
# iteration.df <- data.frame(genomeID = genomeID.table, iterations = itcount.table, 
#                            difference = diffcount.table, stringsAsFactors = FALSE)
# 
# write.table(iteration.df, file = "iterationcount.csv", append = FALSE, sep = ",", 
#             row.names = FALSE, quote = FALSE, col.names = TRUE)

# #######__________________________ comparing ribosomal vs random seed __________________________#######
# 
# 
# # read in both iteration count files, get only first 2 columns
# ribo.seed.iter <- read.csv("itcount_final.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)[, 1:2]
# ribo.seed.iter[ribo.seed.iter == 0] <- NA
# ribo.seed.iter <- na.omit(ribo.seed.iter)
# 
# rand.seed.iter <- read.csv("itcount_rand_final1.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)[, 1:2]
# rand.seed.iter[rand.seed.iter == 0] <- NA
# rand.seed.iter <- na.omit(rand.seed.iter)
# 
# 
# # open files and making suitable for analysis
# # reading a .csv file containing the genome names in the first column
# genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
#                                  as.is=TRUE) #as.is to keep the it as char
# 
# genomeID <- "GCA_000024425"
# 
# # initialize empty colums to add data to it
# genomeID.col <- NULL
# diff.mean.w.col <- NULL
# diff.itcount.col <- NULL
# diff.sum.w.col <- NULL
# for (genomeID in genome.and.organisms[,1]) { 
#   cat (genomeID, "\n")
#   # get all the filenames
#   w.rand.files <- paste("robustnesscheck_iterated_weight_tables/", genomeID, "_robust_it.csv", sep = "")
#   w.ribo.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
#   if (file.exists(w.rand.files)){
#     # open weight files
#     rand.w <- read.csv(file = w.rand.files, sep = ",", header = FALSE, as.is = TRUE)
#     ribo.w <- read.csv(file = w.ribo.files, sep = ",", header = FALSE, as.is = TRUE)
#     # get the difference between codonweights of random seed and ribosomal seed tables
#     diff.w <- rand.w[,2] - ribo.w[,2]
#     # calculate the mean difference
#     diff.mean.w <- mean(diff.w)
#     diff.sum.w <- sum(diff.w)
#     # add the mean difference to the column
#     diff.mean.w.col <- c(diff.mean.w.col, diff.mean.w)
#     diff.sum.w.col <- c(diff.sum.w.col, diff.sum.w)
#     genomeID.col <- c(genomeID.col, genomeID)
#     # get the difference between iteration counts to compute final weight table
#     # when input is ribosomal and random genes
#     if (genomeID %in% rand.seed.iter[,1]){
#       if (genomeID %in% ribo.seed.iter[,1]){
#       itcount.rand <- rand.seed.iter[which(rand.seed.iter[,1] == genomeID), 2]
#       itcount.ribo <- ribo.seed.iter[which(ribo.seed.iter[,1] == genomeID), 2]
#       diff.itcount <- itcount.rand - itcount.ribo
#       diff.itcount.col <- c(diff.itcount.col, diff.itcount)
#       } else { # if the genome IDs are not found fill the column with NA, 
#         # else dataframe will have unequal amount of rows when bound
#         diff.itcount.col <- c(diff.itcount.col, NA)
#         
#       }
#     } else {
#       diff.itcount.col <- c(diff.itcount.col, NA)
#     }
#   }
# }
# 
# data.CAIGCpolII <- data.frame(genomeID.col, diff.mean.w.col, diff.sum.w.col,
#                               diff.itcount.col, stringsAsFactors = FALSE)
# 
# 
# 
# # drawing the histograms 
# # 1st for iteration count difference
# # just the difference
# xlim <- range(data.CAIGCpolII[,3], na.rm = TRUE)
# 
# 
# hist.itcount <- hist(data.CAIGCpolII[,3], breaks=seq(xlim[1], xlim[2]), xlim = xlim, 
#      main = "Difference iteration counts ribosomal vs. random seed", 
#      xlab = "itcount random seed - itcount ribosomal seed", ylab = "number of genomes",
#      col = rgb(0,0,1,1/4))
# 
# # get the amount of genomes per difference
# itcount.diff.hist.info <- data.frame(hist.itcount$breaks[-23], hist.itcount$counts)
# more.it <- sum(itcount.diff.hist.info[1:10,2])/sum(itcount.diff.hist.info[,2])*100
# less.it <- sum(itcount.diff.hist.info[12:22,2])/sum(itcount.diff.hist.info[,2])*100
# no.diff <- itcount.diff.hist.info[11,2]/sum(itcount.diff.hist.info[,2])*100
# 
# xlim <- range(data.CAIGCpolII[,3], na.rm = TRUE)
# 
# 
# hist.itcount <- hist(data.CAIGCpolII[,3], breaks=seq(xlim[1], xlim[2]), xlim = xlim, 
#                      main = "Difference iteration counts ribosomal vs. random seed", 
#                      xlab = "itcount random seed - itcount ribosomal seed", ylab = "number of genomes",
#                      col = rgb(0,0,1,1/4))
# 
# # absolute differnce
# xlim <- range(abs(data.CAIGCpolII[,4]), na.rm = TRUE)
# 
# 
# hist(abs(data.CAIGCpolII[,4]), breaks=seq(xlim[1], xlim[2]), xlim = xlim, 
#      main = "Difference iteration counts ribosomal vs. random seed", 
#      xlab = "abs(itcount random seed - itcount ribosomal seed", ylab = "number of genomes",
#      col = rgb(0,0,1,1/4))
# 
# 
# # drawing the histograms 
# # 2nd for weight table absolute difference
# xlim <- range(abs(data.CAIGCpolII[,2]), na.rm = TRUE)
# nbins <- 8
# 
# 
# hist.wdiff<- hist(abs(data.CAIGCpolII[,2]), breaks=seq(xlim[1], xlim[2], length = nbins), xlim=xlim, 
#      main = "Absolute mean difference codon weights ribosomal vs. random seed", 
#      xlab = "mean(abs(codonweights ribosomal seed - codonweights random seed)", ylab = "number of genomes",
#      col = rgb(0,0,1,1/4))
# hist.wdiff$breaks
# hist.wdiff$count
# 
# 
# 
# 
# # drawing the histograms 
# # 3nd for weight table summed/cumultative difference
# xlim <- range(abs(data.CAIGCpolII[,3]), na.rm = TRUE)
# nbins <- 11
# 
# 
# hist.wdiff<- hist(abs(data.CAIGCpolII[,3]), breaks=seq(xlim[1], xlim[2], length = nbins), xlim=xlim, 
#                   main = "Cumultative difference codon weights ribosomal vs. random seed", 
#                   xlab = "sum(codonweights ribosomal seed - codonweights random seed)", ylab = "number of genomes",
#                   col = rgb(0,0,1,1/4))
# hist.wdiff$breaks
# hist.wdiff$count

