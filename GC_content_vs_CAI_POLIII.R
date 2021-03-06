#!usr/bin/env Rscript

###################################################################
##Author: Melanie van den Bosch
##Script for calculating the mean CAI vs the 
##GC content of all genomes grouped on phyla that are known
##which polIII isoform they contain, then drawing graphs
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
library("ggplot2")
library("grDevices")



#####______________________________________ calculating & writing to file ______________________________________#####

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

#genomeID <- "GCA_000003645"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# read in the polIII data and the bacteria containings these isomers
polIII.data <- read.csv("list_alphasubunitsPOLIII_bacteria.csv", header = TRUE,
                        as.is=TRUE)
polIII.data <- na.omit(polIII.data)


# making the groups which genus is associated to each polIII isomer
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


# filling a dataframe with all the information on meanCAI
# GC content and polIII isomer per genome
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
    # split data on all the polIII groups and add into the dataframe
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
# combine all the columns in a dataframe and write to file
data.CAIGCpolIII <- data.frame(genomeID.col, mean.col, GCcont.col, 
                              polIII.col, stringsAsFactors = FALSE)

fileout = "CAI_GCcont_POLIII_allgenomes.csv"
if (!file.exists(fileout)){
  write.csv(data.CAIGCpolIII, file = fileout, row.names = FALSE)
} else{ 
  print ("file already exists")}
    




#####________________________________draw the plots______________________________#####


data.CAIGCpolIII <- read.table("CAI_GCcont_POLIII_allgenomes.csv", sep = ",", header = TRUE)

# remove genomes that have no polIII isomer
data.CAIGCpolIII[data.CAIGCpolIII[,4] == "Other", 4] <- NA
data.CAIGCpolIII <- na.omit(data.CAIGCpolIII)
    

# draw the plot (normal R)
plot(data.CAIGCpolIII[,3], data.CAIGCpolIII[,2], 
      col = data.CAIGCpolIII$polIII.col,
      main = "Average CAI vs. GC content",
      xlab = "GC content (%)", 
      ylab = "Mean CAI")
grid(NULL, NULL, lty = 6, col = "cornsilk2")
legend("bottomright" ,c("polC/dnaE3", "dnaE1", "dnaE2/dnaE1"), cex=1.5, pch=18,
       col = unique(data.CAIGCpolIII$polIII.col), bty="n")
abline(lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="dnaE1"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="dnaE1"]), col = "black")
abline(lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="dnaE2/dnaE1"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="dnaE2/dnaE1"]), col = "red")
abline(lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="polC/dnaE3"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="polC/dnaE3"]), col = "blue")

# get the data from the regression lines.

lm.dnaE1 <- lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="dnaE1"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="dnaE1"])
eq.dnaE1 <- substitute(y == a + b %.% x*","~~r^2~"="~r2,
                list(a = format(coef(lm.dnaE1)[1], digits = 2),
                     b = format(coef(lm.dnaE1)[2], digits = 2),
                     r2 = format(summary(lm.dnaE1)$r.squared, digits = 3)))

eq.dnaE1 <- as.character(as.expression(eq.dnaE1))

lm.dnaE2.E1 <- lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="dnaE2/dnaE1"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="dnaE2/dnaE1"])
eq.dnaE2 <- substitute(y == a + b %.% x*","~~r^2~"="~r2,
                       list(a = format(coef(lm.dnaE2.E1)[1], digits = 2),
                            b = format(coef(lm.dnaE2.E1)[2], digits = 2),
                            r2 = format(summary(lm.dnaE2.E1)$r.squared, digits = 3)))

eq.dnaE2 <- as.character(as.expression(eq.dnaE2))

lm.polC.dnaE3 <- lm(data.CAIGCpolIII$mean.col[data.CAIGCpolIII$polIII.col=="polC/dnaE3"] ~ 
            data.CAIGCpolIII$GCcont.col[data.CAIGCpolIII$polIII.col=="polC/dnaE3"])
eq.polC <- substitute(y == a + b %.% x*","~~r^2~"="~r2,
                       list(a = format(coef(lm.polC.dnaE3)[1], digits = 2),
                            b = format(coef(lm.polC.dnaE3)[2], digits = 2),
                            r2 = format(summary(lm.polC.dnaE3)$r.squared, digits = 3)))

eq.polC <- as.character(as.expression(eq.polC))


par(mfrow = c(1,1))
# draw the plot (ggplot)
p <- ggplot(data.CAIGCpolIII, aes(y=data.CAIGCpolIII[,2], x= data.CAIGCpolIII[,3], col = data.CAIGCpolIII$polIII.col)) +
  geom_point(size = 2, shape = 18) + 
  geom_smooth(method = "lm", fill = NA)+
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  theme(legend.text = element_text(size = 15))  +
  guides(color = guide_legend(override.aes = list(linetype = c(rep("blank", 3)), 
                                                  size=5))) +
  geom_text(x = 60, y = 0.5, label=eq.dnaE1, parse=TRUE, color = "#F8766D", family = "Helvetica", size = 4.5) + 
  geom_text(x = 60, y = 0.475, label=eq.dnaE2, parse=TRUE, color = "#00BA38", family = "Helvetica", size = 4.5) + 
  geom_text(x = 59, y = 0.45, label=eq.polC, parse=TRUE, color = "#619CFF", family = "Helvetica", size = 4.5) + 
  labs(x = "GC content (%)", y = "Mean CAI",
       color = "Polymerase III subunit isoforms")

p         

# old script
# genomecount = 0
# #pdf("GCvsCAI_plot_polIIIregline.pdf")
#     xlim = c(20, 80)
#     ylim = c(0.35, 0.8)
#     if (genomeID %in% group.dnaE2.dnaE1){
#       col = "blue"
#     } else if (genomeID %in%  group.dnaE1 ){
#         col = "red"
#     } else if(genomeID %in% group.polC.dnaE3) {
#       col = "green"
#     } else {
#       col = rgb(1,1,1,alpha=0.1)
#     }
#     # then plot the GC content against the mean CAI
#       # #if(genomeID %in% biased.genomes){
#       #   print ("biased genome")
#       #   col = "blue"
#       #   } else if (genomeID %in% unbiased.genomes){
#       #     print ("unbiased genome")
#       #   col = "red"
#       #   } else {
#       #   col = "grey"
#       #   }
#     if (genomecount == 0){  
#       plot(GCcont, mean.cai, col=col, 
#                    type = "p", xlim = xlim, ylim = ylim,
#            pch = 18, main = "Average CAI vs. GC content",
#            xlab = "GC content (%)", 
#            ylab = "Mean CAI")
#       grid(NULL, NULL, lty = 6, col = "cornsilk2")
#       legend("bottomright" ,c("dnaE1", "dnaE2/dnaE1", "polC/dnaE3", "Other"), cex=1.5, pch=18,
#             col=c("red", "blue", "green", "grey") , bty="n")
#       genomecount = genomecount + 1
#       } else {
#       points(GCcont, mean.cai, col= col, pch = 18)
#       }
# 
#     }
# }
# dev.off()

#write.table(GCneut, "GCneutralgenomes.csv", sep = ",", 
#            quote = FALSE, col.names = FALSE, row.names = FALSE)
#write.table(GCext, "GCextremegenomes.csv", sep = ",", 
#            quote = FALSE, col.names = FALSE, row.names = FALSE)





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

