#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for visualizing relation between loadings of the PCA, 
## the PC's and extreme GC content genomes
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)
library("Biostrings")
library("RColorBrewer")




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



####________________________ Retrieve data of all the codon weights ________________________####

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

load("GenomeDataSet.RData")


####_______ separating the codons based on aa  _______####

# get the codontable 1 and convert it to a dataframe
aa1 <- getGeneticCode("SGC0")
aatable <- data.frame(keyName = aa1, value = names(aa1), 
                      row.names = NULL, stringsAsFactors = FALSE)
colnames(aatable) <- c("aa", "codon")
# changing aa abbreviations from 1 to 3 letters
aatable[aatable == "A"] <- "Ala"
aatable[aatable == "R"] <- "Arg"
aatable[aatable == "N"] <- "Asn"
aatable[aatable == "D"] <- "Asp"
aatable[aatable == "C"] <- "Cys"
aatable[aatable == "E"] <- "Glu"
aatable[aatable == "Q"] <- "Gln"
aatable[aatable == "G"] <- "Gly"
aatable[aatable == "H"] <- "His"
aatable[aatable == "I"] <- "Ile"
aatable[aatable == "L"] <- "Leu"
aatable[aatable == "K"] <- "Lys"
aatable[aatable == "M"] <- "Met"
aatable[aatable == "F"] <- "Phe"
aatable[aatable == "P"] <- "Pro"
aatable[aatable == "S"] <- "Ser"
aatable[aatable == "T"] <- "Thr"
aatable[aatable == "W"] <- "Trp"
aatable[aatable == "Y"] <- "Tyr"
aatable[aatable == "V"] <- "Val"
aatable[aatable == "*"] <- "Stp"

# and order on aa and then on codon
ordered.aa <- aatable[with(aatable, order(aatable[,1], aatable[,2])), ]

# paste the aa names to the dataframe
codgendf <- data.frame(ordered.aa[,1], codgendf)
#codgendf <- data.frame(rownames(codgendf), codgendf)



# splitting the genomes that are left and right in the PCA plot
length(genomes.rtop.biplot)
length(genomes.ltop.biplot)
codgen.rtop <- codgendf[, genomes.rtop.biplot]
codgen.ltop <- codgendf[, genomes.ltop.biplot]
codgen.rtop <- data.frame(ordered.aa[,1], codgen.rtop)
codgen.ltop <- data.frame(ordered.aa[,1], codgen.ltop)


# setting plotting frame
mar.default <- c(4,3,1,1) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
par(mfrow = c(2,2))
par(cex.axis = 0.5)


# setting custom colors
#palette(rainbow(length(unique(codgendf[,1]))))
#palette(brewer.pal(length(unique(codgendf[,1])), "Spectral")(length(unique(codgendf[,1]))))
n <- length(unique(codgendf[,1]))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
palette(sample(col_vector, n))



# draw boxplots
# right top of PCA
boxplot(t(codgen.rtop[,2:ncol(codgen.rtop)]), col = codgen.rtop[,1],
        xlab = "relative adaptiveness", 
        main = "righttop, should have bias towards ending with C",
        horizontal = TRUE, las = 1)

boxplot(t(codgen.ltop[,2:ncol(codgen.ltop)]), col = codgen.ltop[,1],
        xlab = "relative adaptiveness", 
        main = "lefttop, should have bias towards ending with T/A",
        horizontal = TRUE, las = 1)


# setting plotting frame
mar.default <- c(4,3,1,1) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
par(mfrow = c(2,2))
par(cex.axis = 0.5)

# draw boxplots for genomes based on end position of codon
codgenrtop.G <- codgen.rtop[grepl("G$", rownames(codgen.rtop)), ]
codgenrtop.C <- codgen.rtop[grepl("C$", rownames(codgen.rtop)), ]
codgenrtop.A <- codgen.rtop[grepl("A$", rownames(codgen.rtop)), ]
codgenrtop.T <- codgen.rtop[grepl("T$", rownames(codgen.rtop)), ]

boxplot(t(codgenrtop.C[,2:ncol(codgenrtop.C)]),
        col = "blue",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenrtop.G[,2:ncol(codgenrtop.G)]),
        col = "grey",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenrtop.A[,2:ncol(codgenrtop.A)]),
        col = "green",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenrtop.T[,2:ncol(codgenrtop.T)]),
        col = "red",
        cex.axis = 1.5,
        las = 2)

codgenltop.G <- codgen.ltop[grepl("G$", rownames(codgen.ltop)), ]
codgenltop.C <- codgen.ltop[grepl("C$", rownames(codgen.ltop)), ]
codgenltop.A <- codgen.ltop[grepl("A$", rownames(codgen.ltop)), ]
codgenltop.T <- codgen.ltop[grepl("T$", rownames(codgen.ltop)), ]

boxplot(t(codgenltop.C[,2:ncol(codgenltop.C)]),
        col = "blue",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenltop.G[,2:ncol(codgenltop.G)]),
        col = "grey",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenltop.A[,2:ncol(codgenltop.A)]),
        col = "green",
        cex.axis = 1.5,
        las = 2)
boxplot(t(codgenltop.T[,2:ncol(codgenltop.T)]),
        col = "red",
        cex.axis = 1.5,
        las = 2)
