#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: CAI of organisms in different environmental conditions   
####  Purpose of script: from the gold database different data is retrieved of which
####  environmental conditions a species lives in.
####  The CAI is calculated for these species and compared with each other
#################################################################################

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
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  if (file.exists(cai.files)){
    if (genomeID %in% NeutGen){
      data <- read.csv(file = cai.files, sep = ",", 
                       header = TRUE, as.is = TRUE)
      mean.cai <- mean(data[gencount,2])
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


# draw the plot (normal R)
hist(salinity.df[,3], 
      freq = TRUE,
      col = salinity.df$sal,
      xlab="mean CAI", 
      main="CAI distribution of different salinity requirements", 
      ylab = "Bacterial frequency")
#grid(NULL, NULL, lty = 6, col = "cornsilk2")
legend("topright" ,c("halophile", "halotolerant", "stenohaline", "euryhaline"), cex=1.5, pch=18,
       col = unique(data.CAIGCpolIII$polIII.col), bty="n")

# setting the properties for drawing the graphs
xlim <- range(salinity.df[,3], na.rm = TRUE)
nbins = 50
par(mfrow = c(1,1))

# drawing histogram with all data inside
hist(salinity.df[,3], breaks = seq(xlim[1], xlim[2], length = nbins),
     xlim = xlim)

# retrieving CAI values of each of the conditions
#euryhal.cai <- salinity.df[which(salinity.df$sal == "Euryhaline"), 3]
halophile.cai <- salinity.df[which(salinity.df$sal == "Halophile"), 3]
halotol.cai <- salinity.df[which(salinity.df$sal == "Halotolerant"), 3]
stenohal.cai <- salinity.df[which(salinity.df$sal == "Stenohaline"), 3]

par(mfrow = c(2,2))
# euryhal.hist <- hist(euryhal.cai,
#                      breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

halophile.hist <- hist(halophile.cai,
                      breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)
halophile.hist$density = halophile.hist$counts/sum(halophile.hist$counts)*100

halotol.hist <- hist(halotol.cai,
                       breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)
halotol.hist$density = halotol.hist$counts/sum(halotol.hist$counts)*100


stenohal.hist <- hist(stenohal.cai,
                     breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)
stenohal.hist$density = stenohal.hist$counts/sum(stenohal.hist$counts)*100

par(mfrow = c(2,2))
plot(halophile.hist, freq = FALSE)
plot(halotol.hist, freq = FALSE)
plot(stenohal.hist, freq = FALSE)


# Overlap the histograms

# setting y-axis limit to draw the combined histogram
#ylim <- c(0,max(halophile.hist$density, halotol.hist$density, stenohal.hist$density))
par(mfrow = c(1,1))

plot(halophile.hist, col = rgb(1,1,1/4,1/2), xlim = xlim,  freq=TRUE,
      xlab="mean CAI", cex = 1.5, main="CAI distribution of different salinity requirements",
      ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5)

plot(euryhal.hist, col = rgb(0,0,1,1/4), xlim =  xlim, add = TRUE,  freq = TRUE)
plot(halotol.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = TRUE)
plot(stenohal.hist, col = rgb(0,1,0,1/4), xlim =  xlim, add = TRUE,  freq = TRUE)

legend("topright", c("Halophile", "Euryhaline", "Halotolerant", "Stenohaline"), pch = 15,
       col = c(rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0,1,0,1/4)) , bty = "n", cex=1.5)

# removed eury- and stenohaline
ylim <- c(0,max(halophile.hist$density, halotol.hist$density))


par(mfrow = c(1,1))

plot(halophile.hist, col = rgb(0,0,1,1/4), xlim = xlim,  freq=FALSE,
     xlab="mean CAI", cex = 1.5, main="CAI distribution of different salinity requirements",
     ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5, ylim = ylim)

# plot(euryhal.hist, col = 'yellow', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(halotol.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)
#plot(stenohal.hist, col = 'green', xlim =  xlim, add = TRUE,  freq = FALSE)

legend("topright", c("Halophile", "Halotolerant"), pch = 15,
       col = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), bty = "n", cex=1.5)


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
      mean.cai <- mean(data[gencount,2])
      # setting properties for drawing the graphs
      genID <- c(as.character(genomeID), genID)
      means <- c(mean.cai, means)
      oxy <- c(as.character(oxygen[gencount,2]), oxy)
      gencount = sum(gencount, 1)
    }
  }
}
oxygen.df <- data.frame(genID, oxy, means, stringsAsFactors = FALSE)
length(oxygen.df[,1])



# setting the properties for drawing the graphs
xlim <- range(oxygen.df[,3], na.rm = TRUE)
nbins = 20
par(mfrow = c(1,1))

# drawing histogram with all data inside 
hist(oxygen.df[,3], breaks = seq(xlim[1], xlim[2], length = nbins), 
     xlim = xlim)

# retrieving CAI values of each of the conditions
aerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Aerobe"), 3]
length(aerobe.cai)
anaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Anaerobe"), 3]
length(anaerobe.cai)
facul.cai <- oxygen.df[which(oxygen.df$oxy == "Facultative"), 3]
micr.cai <- oxygen.df[which(oxygen.df$oxy == "Microaerophilic"), 3]
oblaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate aerobe"), 3]
oblanaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate anaerobe"), 3]

par(mfrow = c(2,3))
aerobe.hist <- hist(aerobe.cai, 
                     breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

anearobe.hist <- hist(anaerobe.cai, 
                       breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

facul.hist <- hist(facul.cai, 
                     breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

micr.hist <- hist(micr.cai, 
                      breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

oblaerobe.hist <- hist(oblaerobe.cai, 
                  breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

oblanearobe.hist <- hist(oblanaerobe.cai, 
                  breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)


# Overlap the histograms

# setting y-axis limit to draw the combined histogram
ylim <- c(0.4, max(aerobe.hist$density, anearobe.hist$density, facul.hist$density, 
                micr.hist$density, oblaerobe.hist$density, oblanearobe.hist$density))
par(mfrow = c(1,1))

plot(aerobe.hist, col = rgb(1,1,1/4,1/2), xlim = xlim,  freq=FALSE,
     xlab="mean CAI", cex = 1.5, main="CAI distribution of different oxygen requirements", 
     ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

plot(anearobe.hist, col = rgb(0,0,1,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)
plot(facul.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)
plot(micr.hist, col = rgb(0,1,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)
plot(oblaerobe.hist, col = rgb(0,1,1,1/2), xlim =  xlim, add = TRUE,  freq = FALSE)
plot(oblanearobe.hist, col = rgb(1,1/4,1,1/2), xlim =  xlim, add = TRUE,  freq = FALSE)


legend("topright", c("Aerobe", "Anearobe", "Facultative", "Microareophilic", 
                     "Obligate aerobe", "Obligate anaerobe"), pch = 15,
                      col = c(rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(1,0,0,1/4), 
                              rgb(0,1,0,1/4), rgb(0,1,1,1/2), rgb(1,1/4,1,1/2)), bty = "n", cex=1.5)

# removed eury- and stenohaline
ylim <- c(0,max(halophile.hist$density, halotol.hist$density))


par(mfrow = c(1,1))

plot(halophile.hist, col = rgb(0,0,1,1/4), xlim = xlim,  freq=FALSE,
     xlab="mean CAI", cex = 1.5, main="CAI distribution of different salinity requirements", 
     ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

# plot(euryhal.hist, col = 'yellow', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(halotol.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)
#plot(stenohal.hist, col = 'green', xlim =  xlim, add = TRUE,  freq = FALSE)

legend("topright", c("Halophile", "Halotolerant"), pch = 15,
       col = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), bty = "n", cex=1.5)

