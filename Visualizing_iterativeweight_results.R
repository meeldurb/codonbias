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

# open iterationcount file
itcountdf <- read.csv(file = "itcount_final.csv", header = TRUE, sep=",")
itcountdf <- itcountdf[,1:2]
itcountdf <- na.omit(itcountdf)


# setting the properties for drawing the graph
xlim <- range(itcount[,2], na.rm = TRUE)
nbins = 21


hist(itcount[,2], breaks=seq(xlim[1], xlim[2], length = nbins), xlim = xlim, 
     main = "Iteration counts to compute weight tables for all genomes", 
     xlab = "number of iterations per resulting weight table", ylab = "number of genomes",
     col = rgb(0,0,1,1/4))

#selecting the genomeIDs that are below the cut-off value of 5
# write biased and unbiased genomes to a file
biased.i <- which(itcountdf[,2] < 5)
biased.genomes <- as.vector(itcountdf[biased.i,1])
length(biased.genomes)

unbiased.i <- which(itcountdf[,2] > 5)
unbiased.genomes <- as.vector(itcountdf[unbiased.i,1])
length(unbiased.genomes)

# combine both lists
#un.biased.genomes <- cbind(biased.genomes, unbiased.genomes)

# write to file
write.table(unbiased.genomes, file = "biased.genomes.csv", append=F, 
            sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(unbiased.genomes, file = "unbiased.genomes.csv", append=F, 
            sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)