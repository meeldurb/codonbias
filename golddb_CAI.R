#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: CAI of organisms in different environmental conditions   
####  Purpose of script: from the gold database different data is retrieved of which
####  environmental conditions a species lives in.
####  The CAI is calculated for these species and compared with each other
#################################################################################

install.packages("HistogramTools", repos="http://cran.rstudio.com/")
library("HistogramTools")
library("lattice")


setwd("~/Documents/Master_Thesis_SSB/git_scripts")

genomeID <- "GCA_000003645"
gold.data <- read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
NeutGen <- read.csv(file= "GCneutralgenomes.csv", header = FALSE,
                    as.is = TRUE)
NeutGen <- as.vector(NeutGen[,1])

# what is the information included in this db
names(gold.data)

##############_______________________________Salinity_______________________________#####################
# taking the genomeIDs that have salinity information
salinity <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data)=="Salinity")]))
salinity <- na.omit(salinity)
length(salinity$Salinity)
# check which groups are there and how large they are
unique(salinity$Salinity)
table(salinity[,2])

means = numeric()
genID = character() 
sal = character()
gencount = 1
for (genomeID in salinity[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  if (file.exists(cai.files)){
    if (genomeID %in% NeutGen){
      data <- read.csv(file = cai.files, sep = ",", 
                       header = TRUE, as.is = TRUE)
      mean.cai <- mean(data[,2])
      # setting properties for drawing the graphs
      genID <- c(as.character(genomeID), genID)
      means <- c(mean.cai, means)
      sal <- c(as.character(salinity[gencount,2]), sal)
      gencount = sum(gencount, 1)
    }
  }
}
salinity.df <- data.frame(genID, sal, means, stringsAsFactors = TRUE)
salinity.df <- salinity.df[!salinity.df[,2] == "Undefined", ]


# retrieving CAI values of each of the conditions
euryhal.cai <- salinity.df[which(salinity.df$sal == "Euryhaline"), 3]
halophile.cai <- salinity.df[which(salinity.df$sal == "Halophile"), 3]
halotol.cai <- salinity.df[which(salinity.df$sal == "Halotolerant"), 3]
stenohal.cai <- salinity.df[which(salinity.df$sal == "Stenohaline"), 3]

# draw the histogram plot
xlim <- range(salinity.df[,3], na.rm = TRUE)
breakpoints <- seq(xlim[1], xlim[2], length.out = 20)

# get the histograms and calculate percentages
par(mfrow = c(1,1))
euryhal.hist <- hist(euryhal.cai,
                     breaks = breakpoints, plot = F)
euryhal.hist$counts = euryhal.hist$counts/sum(euryhal.hist$counts)

halophile.hist <- hist(halophile.cai,
                       breaks = breakpoints, plot = F)
halophile.hist$counts = halophile.hist$counts/sum(halophile.hist$counts)

halotol.hist <- hist(halotol.cai,
                     breaks = breakpoints, plot = F)
halotol.hist$counts = halotol.hist$counts/sum(halotol.hist$counts)

stenohal.hist <- hist(stenohal.cai,
                      breaks = breakpoints, plot = F)
stenohal.hist$counts = stenohal.hist$counts/sum(stenohal.hist$counts)

plot(stenohal.hist, col = rgb(0,1,1,1/4), 
     main = "CAI distribution of salinity requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency")

plot(halophile.hist, col = rgb(1,0,0,1/4), add = TRUE)
plot(euryhal.hist, col = rgb(1,1,1/4,1/2), add = TRUE)
plot(halotol.hist, col = rgb(0,0,1,1/4), add = TRUE)
legend("topright", c("Halophile", "Euryhaline", "Halotolerant", "Stenohaline"), pch = 15,
        col = c(rgb(1,0,0,1/4), rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(0,1,1,1/4)) , bty = "n", cex=1.5)



# removed eury- and stenohaline
par(mfrow = c(1,1))
plot(halotol.hist, col = rgb(1,0,0,1/4))
plot(halophile.hist, col = rgb(0,0,1,1/4), add = TRUE)
legend("topleft", c("Halophile", "Halotolerant"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(0,0,1,1/4)) , bty = "n", cex=1.5)

##############_______________________________Oxygen requirement_______________________________#####################
# taking the genomeIDs that have oxygen information
oxygen <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data)=="Oxygen.Requirement")]))
oxygen <- na.omit(oxygen)
length(oxygen$Oxygen.Requirement)

# check which groups are there and how large they are
unique(oxygen$Oxygen.Requirement)
table(oxygen[,2])


means = numeric()
genID = character() 
oxy = character()
gencount = 1
for (genomeID in oxygen[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  if (file.exists(cai.files)){
    if (genomeID %in% NeutGen){
      data <- read.csv(file = cai.files, sep = ",", 
                       header = TRUE, as.is = TRUE)
      mean.cai <- mean(data[,2])
      # setting properties for drawing the graphs
      genID <- c(as.character(genomeID), genID)
      means <- c(mean.cai, means)
      oxy <- c(as.character(oxygen[gencount,2]), oxy)
      gencount = sum(gencount, 1)
    }
  }
}
# write all the data to a dataframe and remove NA rows
oxygen.df <- data.frame(genID, oxy, means, stringsAsFactors = FALSE)
oxygen.df <- na.omit(oxygen.df)
length(oxygen.df[,1])


# retrieving CAI values of each of the conditions
aerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Aerobe"), 3]
length(aerobe.cai)
anaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Anaerobe"), 3]
length(anaerobe.cai)
facul.cai <- oxygen.df[which(oxygen.df$oxy == "Facultative"), 3]
micr.cai <- oxygen.df[which(oxygen.df$oxy == "Microaerophilic"), 3]
length(micr.cai)
oblaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate aerobe"), 3]
oblanaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate anaerobe"), 3]

# draw the histogram plot
xlim <- range(oxygen.df[,3], na.rm = TRUE)
breakpoints <- seq(xlim[1], xlim[2], length.out = 20)

# get the histograms and calculate percentages
par(mfrow = c(1,1))
aerobe.hist <- hist(aerobe.cai,
                     breaks = breakpoints, plot = F)
aerobe.hist$counts = aerobe.hist$counts/sum(aerobe.hist$counts)

anaerobe.hist <- hist(anaerobe.cai,
                       breaks = breakpoints, plot = F)
anaerobe.hist$counts = anaerobe.hist$counts/sum(anaerobe.hist$counts)

facul.hist <- hist(facul.cai,
                     breaks = breakpoints, plot = F)
facul.hist$counts = facul.hist$counts/sum(facul.hist$counts)

micro.hist <- hist(micr.cai,
                      breaks = breakpoints, plot = F)
micro.hist$counts = micro.hist$counts/sum(micro.hist$counts)

oblaerobe.hist <- hist(oblaerobe.cai,
                   breaks = breakpoints, plot = F)
oblaerobe.hist$counts = oblaerobe.hist$counts/sum(oblaerobe.hist$counts)

oblanaerobe.hist <- hist(oblanaerobe.cai,
                   breaks = breakpoints, plot = F)
oblanaerobe.hist$counts = oblanaerobe.hist$counts/sum(oblanaerobe.hist$counts)

plot(facul.hist, col = rgb(1,0,0,1/4),  add = TRUE,
     main = "CAI distribution of oxygen requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency") 
plot(anaerobe.hist, col = rgb(0,0,1,1/4), add = TRUE)
plot(micro.hist, col = rgb(0,1,0,1/4), add = TRUE)
plot(aerobe.hist, col = rgb(1,1,1/4,1/2), add = TRUE)
plot(oblaerobe.hist, col = rgb(0,1,1,1/2), add = TRUE)
plot(oblanaerobe.hist, col = rgb(1,1/4,1,1/2), add = TRUE)


legend("topleft", c("Aerobe", "Anearobe", "Facultative", "Microareophilic", 
                     "Obligate aerobe", "Obligate anaerobe"), pch = 15,
       col = c(rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(1,0,0,1/4), 
               rgb(0,1,0,1/4), rgb(0,1,1,1/2), rgb(1,1/4,1,1/2)), bty = "n", cex=1.5)


plot(stenohal.hist, col = rgb(0,1,1,1/4), 
     main = "CAI distribution of salinity requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency")

plot(halophile.hist, col = rgb(1,0,0,1/4), add = TRUE)
plot(euryhal.hist, col = rgb(1,1,1/4,1/2), add = TRUE)
plot(halotol.hist, col = rgb(0,0,1,1/4), add = TRUE)
legend("topright", c("Halophile", "Euryhaline", "Halotolerant", "Stenohaline"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(0,1,1,1/4)) , bty = "n", cex=1.5)



# removed eury- and stenohaline
par(mfrow = c(1,1))
plot(halotol.hist, col = rgb(1,0,0,1/4))
plot(halophile.hist, col = rgb(0,0,1,1/4), add = TRUE)
legend("topleft", c("Halophile", "Halotolerant"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(0,0,1,1/4)) , bty = "n", cex=1.5)

