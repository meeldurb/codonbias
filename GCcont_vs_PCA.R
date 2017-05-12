#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for visualizing relation between loadings of the PCA, 
## the PC's and extreme GC content genomes
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)



setwd("~/Documents/Master_Thesis_SSB/git_scripts")


load("GenomeDataSet.RData")
colnames(codgendf)
codgen <- t(codgendf)
rownames(codgen)
str(codgen)

gold.data= read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)
#str(gold.data)

# shorten gold.data for genomeIDs we have in codgendf
gold.data <- gold.data[gold.data$GCA %in% rownames(codgen), ]
str(gold.data)
# shorten also codgen data
codgen <- codgen[rownames(codgen) %in% gold.data$GCA, ]
codgen=codgen[gold.data$GCA,]
str(codgen)

# do the PCA
m.codgen <- as.matrix(codgen)
# remove codons that have only value 1, Met and Trp
m.codgen<- m.codgen[, which(colnames(m.codgen)!="ATG")]
m.codgen<- m.codgen[, which(colnames(m.codgen)!="TGG")]

codgen.pca <- prcomp(m.codgen, scale = TRUE)

# draw biplot
par(mar = c(4, 4, 4, 4))
palette(c("White", "Red"))
biplot(codgen.pca, scale = 0)

pca.summary <- summary(codgen.pca)

# get only the genomes that are in top left and top right of the PCA plot
genomes.rtop.biplot <- names(which(pca.summary$x[,1] > 5 & pca.summary$x[,2] > 0))
genomes.ltop.biplot <- names(which(pca.summary$x[,1] < -5 & pca.summary$x[,2] > 0))
#genomes.top.biplot <- names(which(pca.summary$x[,2] > 0))

# read in the table with the GC content of all genomes
data.CAIGCpolIII <- read.table("CAI_GCcont_POLIII_allgenomes.csv", sep = ",", header = TRUE)
GC.content <- data.CAIGCpolIII[,c(1,3)]

rtop.biplot.GC <- as.numeric(GC.content[which(genomes.rtop.biplot %in% GC.content[,1]),2])
ltop.biplot.GC <- as.numeric(GC.content[which(genomes.ltop.biplot %in% GC.content[,1]),2])

mean(rtop.biplot.GC)
mean(ltop.biplot.GC)


