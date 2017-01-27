#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: CAI of organisms in different environmental conditions   
####  Purpose of script: from the gold database different data is retrieved of which
####  environmental conditions a species lives in.
####  The CAI is calculated for these species and compared with each other
#################################################################################

setwd("~/Documents/Master_Thesis_SSB/git_scripts")


gold.data= read.table(file = "gold_gca.tsv", sep="\t" , header=TRUE,
                      row.names=1,as.is=TRUE)

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
for (genomeID in salinity[,1]){
  cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  if (file.exists(cai.files)){
    data <- read.csv(file = cai.files, sep = ",", 
                     header = FALSE, as.is = TRUE)
    mean.cai <- mean(data[,2])
    # setting properties for drawing the graphs
    means <- c(means, mean.cai)
  }
}
salinity.df <- cbind(salinity, means)


# setting the properties for drawing the graphs
xlim <- range(salinity.df[,3], na.rm = TRUE)
nbins = 50
par(mfrow = c(1,1))

# drawing histogram with all data inside 
hist(salinity.df[,3], breaks = seq(xlim[1], xlim[2], length = nbins), 
     xlim = xlim)

# retrieving CAI values of each of the conditions
euryhal.cai <- salinity.df[which(salinity.df$Salinity == "Euryhaline"), 3]
halophile.cai <- salinity.df[which(salinity.df$Salinity == "Halophile"), 3]
halotol.cai <- salinity.df[which(salinity.df$Salinity == "Halotolerant"), 3]
stenohal.cai <- salinity.df[which(salinity.df$Salinity == "Stenohaline"), 3]

par(mfrow = c(2,2))
euryhal.hist <- hist(euryhal.cai, 
                     breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

halophile.hist <- hist(halophile.cai, 
                      breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

halotol.hist <- hist(halotol.cai, 
                       breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

stenohal.hist <- hist(stenohal.cai, 
                     breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)


# Overlap the histograms

# setting y-axis limit to draw the combined histogram
ylim <- c(0,max(halophile.hist$density, halotol.hist$density, stenohal.hist$density))
par(mfrow = c(1,1))

plot(halophile.hist, col = 'blue', xlim = xlim,  freq=FALSE,
      xlab="mean CAI", cex = 1.5, main="CAI distribution of different salinity requirements", 
      ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

plot(euryhal.hist, col = 'yellow', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(halotol.hist, col = 'red', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(stenohal.hist, col = 'green', xlim =  xlim, add = TRUE,  freq = FALSE)

legend("topright", c("Halophile", "Euryhaline", "Halotolerant", "Stenohaline"), pch = 15,
       col = c('blue', 'yellow', 'red', 'green') , bty = "n", cex=1.5)

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
# taking the genomeIDs that have salinity information
oxygen <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data)=="Oxygen.Requirement")]))
oxygen <- na.omit(oxygen)
length(oxygen$Oxygen.Requirement)

# check which groups are there and how large they are
unique(oxygen$Oxygen.Requirement)
table(oxygen[,2])


means = numeric()
for (genomeID in oxygen[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  if (file.exists(cai.files)){
    data <- read.csv(file = cai.files, sep = ",", 
                     header = FALSE, as.is = TRUE)
    mean.cai <- mean(data[,2])
    # setting properties for drawing the graphs
    means <- c(means, mean.cai)
    } 
    # when the file does not exists add a 0 to the means, 
    # else the lenght of the vector will be too small to add to the dataframe
    else{
      means <- c(means, 0)
    }
  }

# the else statement put 0's but we do not want to take these into consideration for our plotting
oxygen.df <- cbind(oxygen, means)
oxygen.df$means[which(oxygen.df$means==0)] <- NA


# setting the properties for drawing the graphs
xlim <- range(oxygen.df[,3], na.rm = TRUE)
nbins = 20
par(mfrow = c(1,1))

# drawing histogram with all data inside 
hist(oxygen.df[,3], breaks = seq(xlim[1], xlim[2], length = nbins), 
     xlim = xlim)

# retrieving CAI values of each of the conditions
aerobe.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Aerobe"), 3]
anaerobe.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Anaerobe"), 3]
facul.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Facultative"), 3]
micr.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Microaerophilic"), 3]
oblaerobe.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Obligate aerobe"), 3]
oblanaerobe.cai <- oxygen.df[which(oxygen.df$Oxygen.Requirement == "Obligate anaerobe"), 3]

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

plot(aerobe.hist, col = 'blue', xlim = xlim,  freq=FALSE,
     xlab="mean CAI", cex = 1.5, main="CAI distribution of different oxygen requirements", 
     ylab = "Bacterial frequency", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

plot(anearobe.hist, col = 'red', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(facul.hist, col = 'green', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(micr.hist, col = 'purple', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(oblaerobe.hist, col = 'orange', xlim =  xlim, add = TRUE,  freq = FALSE)
plot(oblanearobe.hist, col = 'yellow', xlim =  xlim, add = TRUE,  freq = FALSE)


legend("topright", c("Aerobe", "Anearobe", "Facultative", "Microareophilic", 
                     "Obligate aerobe", "Obligate anaerobe"), pch = 15,
                      col = c('blue', 'red', 'green', 'purple', 'orange', 'yellow') , bty = "n", cex=1.5)

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

