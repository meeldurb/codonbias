#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Plotting iteration results of iterative weight tables 
####  Purpose of script: We want to plot the number of iterations for each genome
####  when the weight tables were computed. 
####  depending on the results/pattern seen we will draw conclusions on which genomes
####  are biased or are under translational selection
#################################################################################


setwd("~/Documents/Master_Thesis_SSB/git_scripts")


install.packages(c("ggplot2","RColorBrewer","scales"))
library(ggplot2); library(scales); library(grid); library(RColorBrewer)


# open iterationcount file
itcountdf <- read.csv(file = "itcount_final.csv", header = TRUE, sep=",")
itcountdf <- itcountdf[,1:2]
itcountdf <- na.omit(itcountdf)


# setting the properties for drawing the graph
xlim <- range(itcountdf[,2], na.rm = TRUE)
nbins = 21


ggplot(data = itcountdf, aes(itcountdf[,2])) +
      geom_histogram(breaks = seq(0, 21, by=0.5),
                     col = "black") +
      labs(x = "number of iterations per resulting weight table", 
           y = "number of genomes")

hist <- hist(itcountdf[,2], breaks=seq(0, 21, 1), xlim = xlim, 
     main = NA, 
     xlab = "number of iterations per resulting weight table", ylab = "number of genomes",
     col = rgb(0,0,1,1/4))

#selecting the genomeIDs that are below the cut-off value of 5
biased.i <- which(itcountdf[,2] <= 4 & itcountdf[,2] != 0)
biased.genomes <- as.vector(itcountdf[biased.i,1])

unbiased.i <- which(itcountdf[,2] >= 5)
unbiased.genomes <- as.vector(itcountdf[unbiased.i,1])

# numbers of genomes contained in the datasets
length(which(itcountdf[,2]==0))
length(biased.genomes)
length(unbiased.genomes)

length(which(itcountdf[,2]==0)) + length(biased.genomes) + length(unbiased.genomes) 
#should be 6041

# combine both lists
#un.biased.genomes <- cbind(biased.genomes, unbiased.genomes)

# write to file
write.table(biased.genomes, file = "it1_4genomes.csv", append=F, 
            sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(unbiased.genomes, file = "it5_21genomes.csv", append=F, 
            sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
