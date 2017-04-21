#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Visualizing the differences between ribosomal seed and
####  random seed of the iterative weight table algorithm
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



#######__________________________ comparing ribosomal vs random seed __________________________#######


# read in both iteration count files, get only first 2 columns
ribo.seed.iter <- read.csv("itcount_final.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)[, 1:2]
ribo.seed.iter[ribo.seed.iter == 0] <- NA
ribo.seed.iter <- na.omit(ribo.seed.iter)

rand.seed.iter <- read.csv("itcount_rand_final1.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)[, 1:2]
rand.seed.iter[rand.seed.iter == 0] <- NA
rand.seed.iter <- na.omit(rand.seed.iter)


# open files and making suitable for analysis
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE,
                                 as.is=TRUE) #as.is to keep the it as char

genomeID <- "GCA_000024425"

# initialize empty colums to add data to it
genomeID.col <- NULL
diff.mean.w.col <- NULL
diff.itcount.col <- NULL
diff.sum.w.col <- NULL
for (genomeID in genome.and.organisms[,1]) {
  cat (genomeID, "\n")
  # get all the filenames
  w.rand.files <- paste("robustnesscheck_iterated_weight_tables/", genomeID, "_robust_it.csv", sep = "")
  w.ribo.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
  if (file.exists(w.rand.files)){
    # open weight files
    rand.w <- read.csv(file = w.rand.files, sep = ",", header = FALSE, as.is = TRUE)
    ribo.w <- read.csv(file = w.ribo.files, sep = ",", header = FALSE, as.is = TRUE)
    # get the difference between codonweights of random seed and ribosomal seed tables
    diff.w <- rand.w[,2] - ribo.w[,2]
    # calculate the mean difference
    diff.mean.w <- mean(diff.w)
    diff.sum.w <- sum(diff.w)
    # add the mean difference to the column
    diff.mean.w.col <- c(diff.mean.w.col, diff.mean.w)
    diff.sum.w.col <- c(diff.sum.w.col, diff.sum.w)
    genomeID.col <- c(genomeID.col, genomeID)
    # get the difference between iteration counts to compute final weight table
    # when input is ribosomal and random genes
    if (genomeID %in% rand.seed.iter[,1]){
      if (genomeID %in% ribo.seed.iter[,1]){
        itcount.rand <- rand.seed.iter[which(rand.seed.iter[,1] == genomeID), 2]
        itcount.ribo <- ribo.seed.iter[which(ribo.seed.iter[,1] == genomeID), 2]
        diff.itcount <- itcount.rand - itcount.ribo
        diff.itcount.col <- c(diff.itcount.col, diff.itcount)
      } else { # if the genome IDs are not found fill the column with NA,
        # else dataframe will have unequal amount of rows when bound
        diff.itcount.col <- c(diff.itcount.col, NA)
        
      }
    } else {
      diff.itcount.col <- c(diff.itcount.col, NA)
    }
  }
}

data.itcount <- data.frame(genomeID.col, diff.mean.w.col, diff.sum.w.col,
                              diff.itcount.col, stringsAsFactors = FALSE)



# drawing the histograms
# 1st for iteration count difference
# just the difference
xlim <- range(data.itcount[,4], na.rm = TRUE)


hist.itcount <- hist(data.itcount[,4], breaks = seq(xlim[1], xlim[2], by=1),
                     main = "Difference iteration counts ribosomal vs. random seed",
                     xlab = "itcount random seed - itcount ribosomal seed", ylab = "number of genomes",
                     col = rgb(0,0,1,1/4))


ggplot(data = data.itcount, aes(data.itcount[,4])) +
  geom_histogram(breaks = seq(xlim[1], xlim[2], by=1),
                 col = "black") +
  labs(x = "itcount random seed - itcount ribosomal seed", 
       y = "number of genomes")


# get the amount of genomes per difference
itcount.diff.hist.info <- data.frame(hist.itcount$breaks[-18], hist.itcount$counts)
more.it <- sum(itcount.diff.hist.info[1:11,2])/sum(itcount.diff.hist.info[,2])*100
less.it <- sum(itcount.diff.hist.info[13:17,2])/sum(itcount.diff.hist.info[,2])*100
no.diff <- itcount.diff.hist.info[12,2]/sum(itcount.diff.hist.info[,2])*100
sum(itcount.diff.hist.info[,2])




# absolute differnce
xlim <- range(abs(data.itcount[,4]), na.rm = TRUE)


hist(abs(data.itcount[,4]), breaks=seq(xlim[1], xlim[2]), xlim = xlim,
     main = "Difference iteration counts ribosomal vs. random seed",
     xlab = "abs(itcount random seed - itcount ribosomal seed", ylab = "number of genomes",
     col = rgb(0,0,1,1/4))

ggplot(data = data.itcount, aes(abs(data.itcount[,4]))) +
  geom_histogram(breaks = seq(0, 19, by=1),
                col = "black") +
  labs(x = "abs(itcount random seed - itcount ribosomal seed)", 
       y = "number of genomes")

# drawing the histograms
# 2nd for weight table absolute difference
xlim <- range(abs(data.itcount[,2]), na.rm = TRUE)
nbins <- 8


hist.wdiff<- hist(abs(data.itcount[,2]), breaks=seq(xlim[1], xlim[2], length = nbins), xlim=xlim,
                  main = "Absolute mean difference codon weights ribosomal vs. random seed",
                  xlab = "mean(abs(codonweights ribosomal seed - codonweights random seed)", ylab = "number of genomes",
                  col = rgb(0,0,1,1/4))
hist.wdiff$breaks
sum(hist.wdiff$count)




# drawing the histograms
# 3nd for weight table summed/cumultative difference
xlim <- range(abs(data.itcount[,3]), na.rm = TRUE)
nbins <- 11


hist.wdiff<- hist(abs(data.itcount[,3]), breaks=seq(xlim[1], xlim[2], length = nbins), xlim=xlim,
                  main = "Cumultative difference codon weights ribosomal vs. random seed",
                  xlab = "sum(codonweights ribosomal seed - codonweights random seed)", ylab = "number of genomes",
                  col = rgb(0,0,1,1/4))
hist.wdiff$breaks
hist.wdiff$count

