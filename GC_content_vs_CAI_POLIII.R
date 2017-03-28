#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the mean CAI vs the 
##GC content of all genomes grouped on phyla that are known
##which polIII isoform they contain
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



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

#genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# the 3 different groups of polIII isoforms and its groups belonging to it

polIII.data <- read.csv("list_alphasubunitsPOLIII_bacteria.csv", header = TRUE,
                        as.is=TRUE)
polIII.data <- na.omit(polIII.data)


genus.dnaE1 <- NULL
genus.polC.dna3 <- NULL
genus.dnaE1.dnaE2 <- NULL
for (row in seq_along(polIII.data[,1])){
  print (polIII.data[row,16])
  if (polIII.data[row,2] == "" && polIII.data[row,3] == "" && polIII.data[row,4] == "") {
    genus.dnaE1 <- c(genus.dnaE1, polIII.data[row,16]) 
    } 
  else if (polIII.data[row,1] == "" && polIII.data[row,2] == "") {
    genus.polC.dna3 <- c(genus.polC.dna3, polIII.data[row,16])
  }
  else if (polIII.data[row,3] == "" && polIII.data[row,4] == "" && polIII.data[row,2] == "dnaE_II") {
    genus.dnaE1.dnaE2 <- c(genus.dnaE1.dnaE2, polIII.data[row,16])
  }
} 
  
dnaE1table <- table(genus.dnaE1)
dnaE1E2 <- table(genus.dnaE1.dnaE2)
dnapolcE2 <- table(genus.polC.dna3)


uniq.genus.dnaE1 <- unique(genus.dnaE1)
uniq.genus.polC.dna3 <- unique(genus.polC.dna3)
uniq.genus.dnaE1.dnaE2 <- unique(genus.dnaE1.dnaE2)


which(genus.dnaE1 %in% genus.dnaE1.dnaE2)
which(genus.dnaE1.dnaE2 %in% genus.dnaE1)

# duplicated gena?
genus.dnaE1 %in% genus.dnaE1.dnaE2
genus.dnaE1 %in% genus.polC.dna3
genus.dnaE1.dnaE2 %in% genus.polC.dna3

genus.dnaE1
genus.polC.dna3
genus.dnaE1.dnaE2



# Taking the information on phyla from the gold db
gold.data= read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)

# taking the genomeIDs that have genus information
genus <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data)=="NCBI.Genus")]))
genus <- na.omit(genus)
length(genus)
# check which groups are there and how large they are
unique(genus$NCBI.Genus)
table(genus[,2])

# making a list of which of the organisms contain a certain polIII isoform
group.dnaE1 <- genus[which(genus$NCBI.Genus %in% uniq.genus.dnaE1), 1]
length(group.dnaE1)
group.dnaE2.dnaE1 <- genus[which(genus$NCBI.Genus %in% uniq.genus.dnaE1.dnaE2), 1]
length(group.dnaE2.dnaE1)
group.polC.dnaE3 <- genus[which(genus$NCBI.Genus %in% uniq.genus.polC.dna3), 1]
length(group.polC.dnaE3)

length(group.dnaE1) + length(group.dnaE2.dnaE1) + length(group.polC.dnaE3)

# do the lm fitting

# dnaE1
lm.fit.mt = lm(group.dnaE1 ~ )
summary(lm.fit.mt)
par(mfrow = c(2,2))
plot(lm.fit.mt)
mean(lm.fit.mt$residuals^2)



#GCneut <- NULL
#GCext <- NULL
genomeID.col <- NULL
mean.col <- NULL
GCcont.col <- NULL
polIII.col <- NULL

for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  CDS.files <- paste("CDS_data/", genomeID, "_CDS.csv", sep = "")
  if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, sep = ",", header = TRUE, as.is = TRUE)
    seq.data <- read.csv(file = CDS.files, sep = ",", header = TRUE, as.is = TRUE)
    
    genomeID.col <- c(genomeID.col, genomeID)
    # calculate mean of all the CAI values in the genome
    mean.cai <- mean(cai.data[,2])
    mean.col <- c(mean.col, as.numeric(mean.cai))
    # calculate the GC content of the genome
    all.seq <- paste(as.matrix(seq.data)[,2], sep="", collapse="")
    seq.split <- strsplit(all.seq, "")[[1]]
    GCcont <- GC(seq.split)*100
    GCcont.col <- c(GCcont.col, as.numeric(GCcont))
    if (genomeID %in% group.dnaE2.dnaE1){
      polIII.col <- c(polIII.col, "dnaE2/dnaE1")
    } else if (genomeID %in%  group.dnaE1 ){
      polIII.col <- c(polIII.col, "dnaE1")
    } else if(genomeID %in% group.polC.dnaE3) {
      polIII.col <- c(polIII.col, "polC/dnaE3")
    } else {
      polIII.col <- c(polIII.col, "Other")
    }
  }
}
data.CAIGCpolII <- data.frame(genomeID.col, mean.col, GCcont.col, 
                              polIII.col, stringsAsFactors = FALSE)

write.csv(data.CAIGCpolII, file = "CAI_GCcont_POLIII_allgenomes.csv", row.names = FALSE)
    

    
    
    
genomecount = 0
pdf("GCvsCAI_plot_polIIIregline.pdf")
    # if (GCcont > 65.0  | GCcont < 35.0){
    #   GCext <- c(GCext, genomeID)
    # }
    # if (GCcont > 35.0 && GCcont < 65.0){
    #   GCneut <- c(GCneut, genomeID)
    # }
    xlim = c(20, 80)
    ylim = c(0.35, 0.8)
    if (genomeID %in% group.dnaE2.dnaE1){
      col = "blue"
    } else if (genomeID %in%  group.dnaE1 ){
        col = "red"
    } else if(genomeID %in% group.polC.dnaE3) {
      col = "green"
    } else {
      col = rgb(1,1,1,alpha=0.1)
    }
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
      plot(GCcont, mean.cai, col=col, 
                   type = "p", xlim = xlim, ylim = ylim,
           pch = 18, main = "Average CAI vs. GC content",
           xlab = "GC content (%)", 
           ylab = "Mean CAI")
      grid(NULL, NULL, lty = 6, col = "cornsilk2")
      legend("bottomright" ,c("dnaE1", "dnaE2/dnaE1", "polC/dnaE3", "Other"), cex=1.5, pch=18,
            col=c("red", "blue", "green", "grey") , bty="n")
      genomecount = genomecount + 1
      } else {
      points(GCcont, mean.cai, col= col, pch = 18)
      }

    }
}
dev.off()

write.table(GCneut, "GCneutralgenomes.csv", sep = ",", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(GCext, "GCextremegenomes.csv", sep = ",", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)





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

