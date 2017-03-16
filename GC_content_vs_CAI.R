#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the mean CAI vs the 
##GC content of all genomes
###################################################################

# Packages needed to be installed to calculate the frequency of oligonucleotides
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
# to compute CAI
install.packages("seqinr", repos="http://cran.rstudio.com/")
install.packages("locfit", repos="http://cran.rstudio.com/")


# loading required library to compute the codon frequency
library("Biostrings")
# loading required library to use CAI function
library("seqinr")
library("locfit")
# install packages to draw plots
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# open iterationcount file
itcountdf <- read.csv(file = "itcount_final.csv", header = TRUE, sep=",")
itcountdf <- itcountdf[,1:2]
itcountdf <- na.omit(itcountdf)


# getting groups of it counts and factorizing to later colour
sort.count <- rle(sort(itcountdf[,2]))

#order.count <- rle(sort(gold.data$NCBI.Order))
selected <- sort.count$values[which(sort.count$length>1)]   
##put colors in those for which we have more than 5 (increase for a "real" example)
group <- itcountdf$difference.count
group[which(!group%in% selected)] <- NA
itcountdf$Group <- group
itcountdf$Group <- factor(itcountdf$Group, levels=c(selected))


#itcountdf$Group <- factor(itcountdf[,2], levels = itcountdf[,2])
#title <- "ggplot of Orders"

myplot <- ggplot(df, aes(PC1.codgen, PC2.codgen, color=Group))+   #these commands creat the plot, but nothing appears
  geom_point(size=2, shape=18)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 0,  size = 12, hjust = 0.5))+
  theme_bw()+
  theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) +
  theme(legend.text = element_text(size = 10))  +
  xlab(paste("PC1 (", format(pca.summary$importance[2,1]*100, digits = 2),"%)", sep = "")) +   
  ylab(paste("PC2 (", format(pca.summary$importance[2,2]*100, digits = 2),"%)", sep = "")) +
  ggtitle(title)

genomecount = 0
#pdf("GCvsCAI_plot_itcount.pdf")
for (genomeID in itcountdf[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  CDS.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
  if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, sep = ",", header = TRUE, as.is = TRUE)
    seq.data <- read.csv(file = CDS.files, sep = ",", header = TRUE, as.is = TRUE)
    
    # calculate mean of all the CAI values in the genome
    mean.cai <- mean(cai.data[,2])
    
    # calculate the GC content of the genome
    all.seq <- paste(as.matrix(seq.data)[,2], sep="", collapse="")
    seq.split <- strsplit(all.seq, "")[[1]]
    GCcont <- GC(seq.split)*100
    xlim = c(20, 80)
    ylim = c(0.35, 0.8)
    #col = "blue"
    # then plot the GC content against the mean CAI
      # #if(genomeID %in% biased.genomes){
      #   print ("biased genome")
      #   col = "blue"
      #   } else if (genomeID %in% unbiased.genomes){
      #     print ("unbiased genome")
      #   col = "red"
      #   } else {
      #   col = "grey"
      #   }
    if (genomecount == 0){  
      plot(GCcont, mean.cai, col=itcountdf$Group, 
                   type = "p", xlim = xlim, ylim = ylim,
           pch = 18, main = "Average CAI vs. GC content",
           xlab = "GC content (%)", 
           ylab = "Mean CAI")
      grid(NULL, NULL, lty = 6, col = "cornsilk2")
      #legend("topleft" ,c("biased", "unbiased", "unknown"), cex=1.5, pch=18,
            # col=c("blue", "red", "grey") , bty="n")
      genomecount = genomecount + 1
      } else {
      points(GCcont, mean.cai, col= itcountdf$Group, pch = 18)
      }
  }
}

dev.off()



######______________________________For 1 genome______________________________######

genomeID <- "GCA_000003645"
cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
CDS.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
cai.data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)
seq.data <- read.csv(file = CDS.files, sep = ",", header = TRUE, as.is = TRUE)

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


mean.cai <- mean(cai.data[,2])


all.seq <- paste(as.matrix(seq.data)[,2], sep="", collapse="")
GC(all.seq)
seq.split <- strsplit(all.seq, "")[[1]]
GCcont <- GC(seq.split)*100
#nfreq <- table(seq.split)
#GCcont <- ((nfreq[2] + nfreq[3])/sum(nfreq))*100

