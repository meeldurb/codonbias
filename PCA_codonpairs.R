#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script calculating and drawing PCA plots from the codonpair 
##weight tables of all the genomes
###################################################################


# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# open files and making suitable for analysis
genomeID <- "GCA_000008185"
genomeID <- "GCA_000008205"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                as.is=TRUE) #as.is to keep the it as char
genome.and.organisms <- genome.and.organisms[,1:2]



####________________________ Do the PCA of all genomes and codon weights ________________________####

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

load("codonpairsGenomeDataSet.RData")
colnames(codgendf)
codgen <- t(codgendf)
rownames(codgen)
str(codgen)
codgen[is.na(codgen)] <- 0

gold.data= read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)
str(gold.data)

# shorten gold.data for genomeIDs we have in codgendf
gold.data <- gold.data[gold.data$GCA %in% rownames(codgen), ]
str(gold.data)
# shorten also codgen data
codgen <- codgen[rownames(codgen) %in% gold.data$GCA, ]
codgen=codgen[gold.data$GCA,]
str(codgen)

# do the PCA
m.codgen <- as.matrix(codgen)
#m.codgen<- m.codgen[, which(colnames(m.codgen)!="ATG")]

codgen.pca <- prcomp(m.codgen, scale = TRUE)

par(mar = c(5, 4, 4, 4))
palette(c("White", "Red"))
biplot(codgen.pca, scale = 0)

pca.summary <- summary(codgen.pca)

#easy plot
PC1.codgen <- as.numeric(codgen.pca$x[,1])
PC2.codgen <- as.numeric(codgen.pca$x[,2])

plot(PC1.codgen, PC2.codgen, xlab=paste("PC1 (", format(pca.summary$importance[2,1]*100, digits=2),"%)", sep=""), 
     ylab=paste("PC2 (", format(pca.summary$importance[2,2]*100, digits=2),"%)", sep=""))




### ggplot Genus
df <-data.frame(PC1.codgen, PC2.codgen)

order.count <- rle(sort(gold.data$NCBI.Genus))
selected <- order.count$values[which(order.count$length>15)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$NCBI.Genus
group[which(!group%in% selected)] <- "Other"

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
legendtitle <- "Genera"


### ggplot FAMILY
df <-data.frame(PC1.codgen, PC2.codgen)

order.count <- rle(sort(gold.data$NCBI.Family))
selected <- order.count$values[which(order.count$length>20)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$NCBI.Family
group[which(!group%in% selected)] <- "Other"

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
legendtitle <- "Families"




### ggplot ORDER
df <-data.frame(PC1.codgen, PC2.codgen)

order.count <- rle(sort(gold.data$NCBI.Order))
selected <- order.count$values[which(order.count$length>25)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$NCBI.Order
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
legendtitle <- "Orders"

### ggplot CLASS
df <-data.frame(PC1.codgen, PC2.codgen)

class.count <- rle(sort(gold.data$NCBI.Class))
selected <- class.count$values[which(class.count$length>15)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$NCBI.Class
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
legendtitle <- "Classes"

### ggplot PHYLUM
df <-data.frame(PC1.codgen, PC2.codgen)

phylum.count <- rle(sort(gold.data$NCBI.Phylum))
selected <- phylum.count$values[which(phylum.count$length>10)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$NCBI.Phylum
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected) )

legendtitle <- "Phyla"


### ggplot SHAPE
df <-data.frame(PC1.codgen, PC2.codgen)

shape.count <- rle(sort(gold.data$Cell.Shape))
selected <- shape.count$values[which(shape.count$length>1)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$Cell.Shape
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected) )
title <- "Principal components analysis of codonpairs weight in organim shape"

### ggplot GRAM STAIN
df <-data.frame(PC1.codgen, PC2.codgen)

gram.count <- rle(sort(gold.data$Gram.Stain))
selected <- gram.count$values[which(gram.count$length>1)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$Gram.Stain
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
title <- "Principal components analysis of codonpairs weight in gram stain"

### ggplot OXYGEN REQUIREMENT
df <-data.frame(PC1.codgen, PC2.codgen)

oxygen.count <- rle(sort(gold.data$Oxygen.Requirement))
selected <- oxygen.count$values[which(oxygen.count$length>1)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- gold.data$Oxygen.Requirement
group[which(!group%in% selected)] <- NA

df$Group <- group
df$Group <- factor(df$Group, levels=c(selected))
title <- "Principal components analysis of codonpairs weight in oxygen requirement"


### Draw plot
myplot <- ggplot(df, aes(PC1.codgen, PC2.codgen, color=Group))+   #these commands create the plot, but nothing appears
  geom_point(size=4, shape=18)+
  theme_bw(base_size = 13)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5))+
  theme_bw()+
  theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) +
  theme(legend.text = element_text(size = 13))  +
  xlab(paste("PC1 (", format(pca.summary$importance[2,1]*100, digits = 2),"%)", sep = "")) +   
  ylab(paste("PC2 (", format(pca.summary$importance[2,2]*100, digits = 2),"%)", sep = "")) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  scale_colour_discrete(name = legendtitle) +
  scale_fill_brewer(palette = "Spectral")
  #ggtitle(title) 

print(myplot)   #have a look at the plot 



ggsave(file="test.pdf", myplot, width=8, height=8)  #save to a file (change extensdion for tiff, png..)



#another aspect  #DEFAULS
myplot <- ggplot(df,aes( PC1,PC2, color=Group))  +    #these commands creat the plot, but nothing appears
  geom_point(size=3, shape=15)

print(myplot) 

###_______________________________example script from Maria_______________________________###
setwd("~/Documents/Master_Thesis_SSB/git_scripts/old or unused scripts and data/PCA example")

load("toyDataSet.RData")   ##this loads the data 
str(data)   ###100 samples (genomes) and 70 variables (domains)

load("toyMetaDataSet.RData")  ##this loads the metadata
str(domdata) #metadata for the genomes



##do the PCA
yy <- as.matrix(data)
xx <- prcomp(yy)
##get extra info 
z=summary(xx)

#easy plot
PC1 <- as.numeric(xx$x[,1])
PC2 <- as.numeric(xx$x[,2])

plot(PC1, PC2, xlab=paste("PC1 (", format(z$importance[2,1]*100, digits=2),"%)", sep=""), 
     ylab=paste("PC2 (", format(z$importance[2,2]*100, digits=2),"%)", sep=""))


###Nicer   Plot
df <-data.frame(PC1, PC2)

kk <- rle(sort(domdata$Genus ))
selected <- kk$values[which(kk$lengths>5)]   ##put colors in those for which we have more than 5 (increase for a "real" example
group <- domdata$Genus
group[which(!group%in% selected)] <- "Other"

library(ggplot2)
df$Group <-    group
df$Group <- factor(df$Group, levels=c(selected, "Other") )


myplot <- ggplot(df,aes( PC1,PC2, color=Group))+   #these commands creat the plot, but nothing appears
    geom_point(size=3, shape=15)+
    theme(axis.text.x = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5)) + 
    theme(axis.text.y = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5))+
    theme_bw()+
    theme(legend.background = element_rect(fill="white", size=.0, linetype="dotted")) +
    theme(legend.text = element_text( size = 10))  +
    xlab( paste("PC1 (", format(z$importance[2,1]*100, digits=2),"%)", sep="")) +   
    ylab( paste("PC2 (", format(z$importance[2,2]*100, digits=2),"%)", sep=""))



print(myplot)   #have a look at the plot 

ggsave(file="test.pdf", myplot, width=8, height=8)  #save to a file (change extensdion for tiff, png..)



#another aspect  #DEFAULS
myplot <- ggplot(df,aes( PC1,PC2, color=Group))  +    #these commands creat the plot, but nothing appears
    geom_point(size=3, shape=15)

print(myplot) 
