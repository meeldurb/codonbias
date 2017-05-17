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
     cex.axis = 1,
     ylab = "Bacterial frequency")

plot(halophile.hist, col = rgb(1,0,0,1/4), add = TRUE)
plot(euryhal.hist, col = rgb(1,1,1/4,1/2), add = TRUE)
plot(halotol.hist, col = rgb(0,0,1,1/4), add = TRUE)
legend("topright", c("Halophile", "Euryhaline", "Halotolerant", "Stenohaline"), pch = 15,
        col = c(rgb(1,0,0,1/4), rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(0,1,1,1/4)) , bty = "n", cex=1.5)



# removed eury- and stenohaline
par(mfrow = c(1,1))
plot(halotol.hist, col = rgb(1,0,0,1/4), 
     main = "CAI distribution of salinity requirements",
     xlab = "mean CAI", 
     ylab = "Bacterial frequency",
     cex.axis = 1)
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
length(facul.cai)
micr.cai <- oxygen.df[which(oxygen.df$oxy == "Microaerophilic"), 3]
length(micr.cai)
oblaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate aerobe"), 3]
length(oblaerobe.cai)
oblanaerobe.cai <- oxygen.df[which(oxygen.df$oxy == "Obligate anaerobe"), 3]
length(oblanaerobe.cai)


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

plot(facul.hist, col = rgb(1,0,0,1/4),
     main = "CAI distribution of oxygen requirement",
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



# splitting oxygen requirement in 2 groups, one that can live without one that cannot live without
noox <- c(oblanaerobe.cai, facul.cai, anaerobe.cai)
proox <- c(micr.cai, aerobe.cai, oblaerobe.cai)


# draw the histogram plot
xlim <- range(oxygen.df[,3], na.rm = TRUE)
breakpoints <- seq(xlim[1], xlim[2], length.out = 20)

# get the histograms and calculate percentages
par(mfrow = c(1,1))
noox.hist <- hist(noox,
                    breaks = breakpoints, plot = F)
noox.hist$counts = noox.hist$counts/sum(noox.hist$counts)

proox.hist <- hist(proox,
                  breaks = breakpoints, plot = F)
proox.hist$counts = proox.hist$counts/sum(proox.hist$counts)


plot(proox.hist, col = rgb(1,0,0,1/4),
     main = "CAI distribution of oxygen requirement",
     xlab = "mean CAI",
     ylab = "Bacterial frequency",
     cex.axis = 1) 
plot(noox.hist, col = rgb(0,0,1,1/4), add = TRUE)

legend("topleft", c("no oxygen", "pro oxygen"), pch = 15,
       col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = "n", cex=1.5)



##############_______________________________temperature requirement_______________________________#####################
# taking the genomeIDs that have oxygen information
temp <- as.data.frame(c(gold.data[1], gold.data[which(colnames(gold.data) == "Temperature.Range")]))
temp <- na.omit(temp)
length(temp$Temperature.Range)
# check which groups are there and how large they are
unique(temp$Temperature.Range)
table(temp[,2])

# changing the temperature classifiers because of no uniformity
also.meso <- c("25-37", "26", "28", "28 - 30 C", "28-32", "28 C", 
              "30", "30 C", "35 C", "Apr-37", "May-37", "mesophile")
temp$Temperature.Range <- as.character(temp$Temperature.Range)
for (elem in also.meso){
  temp$Temperature.Range[temp$Temperature.Range == elem] <- "Mesophile"
}
temp$Temperature.Range[temp$Temperature.Range == "15-20 C, Psychrophile"] <- "Psychrophile"
temp$Temperature.Range[temp$Temperature.Range == "to 93 C"] <- "Thermophile"

means = numeric()
genID = character() 
t = character()
gencount = 1
for (genomeID in temp[,1]){
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
      t <- c(as.character(temp[gencount,2]), t)
      gencount = sum(gencount, 1)
    }
  }
}
# write all the data to a dataframe and remove NA rows
temp.df <- data.frame(genID, t, means, stringsAsFactors = FALSE)
temp.df <- na.omit(temp.df)
length(temp.df[,1])

# retrieving CAI values of each of the conditions
meso.cai <- temp.df[which(temp.df$t == "Mesophile"), 3]
length(meso.cai)
hypther.cai <- temp.df[which(temp.df$t == "Hyperthermophile"), 3]
length(hypther.cai)
psychro.cai <- temp.df[which(temp.df$t == "Psychrophile"), 3]
length(psychro.cai)
psychrotol.cai <- temp.df[which(temp.df$t == "Psychrotolerant"), 3]
length(psychrotol.cai)
psychrotrop.cai <- temp.df[which(temp.df$t == "Psychrotrophic"), 3]
length(psychrotrop.cai)
thermo.cai <- temp.df[which(temp.df$t == "Thermophile"), 3]
length(thermo.cai)
thermotol.cai <- temp.df[which(temp.df$t == "Thermotolerant"), 3]
length(thermotol.cai)


# draw the histogram plot
xlim <- range(temp.df[,3], na.rm = TRUE)
breakpoints <- seq(xlim[1], xlim[2], length.out = 20)

# get the histograms and calculate percentages
par(mfrow = c(1,1))
meso.hist <- hist(meso.cai,
                    breaks = breakpoints, plot = F)
meso.hist$counts = meso.hist$counts/sum(meso.hist$counts)

hypther.hist <- hist(hypther.cai,
                      breaks = breakpoints, plot = F)
hypther.hist$counts = hypther.hist$counts/sum(hypther.hist$counts)

psychro.hist <- hist(psychro.cai,
                   breaks = breakpoints, plot = F)
psychro.hist$counts = psychro.hist$counts/sum(psychro.hist$counts)

psychrotol.hist <- hist(psychrotol.cai,
                   breaks = breakpoints, plot = F)
psychrotol.hist$counts = psychrotol.hist$counts/sum(psychrotol.hist$counts)

psychrotrop.hist <- hist(psychrotrop.cai,
                       breaks = breakpoints, plot = F)
psychrotrop.hist$counts = psychrotrop.hist$counts/sum(psychrotrop.hist$counts)

thermo.hist <- hist(thermo.cai,
                         breaks = breakpoints, plot = F)
thermo.hist$counts = thermo.hist$counts/sum(thermo.hist$counts)

thermotol.hist <- hist(thermotol.cai,
                    breaks = breakpoints, plot = F)
thermotol.hist$counts = thermo.hist$counts/sum(thermotol.hist$counts)



plot(thermo.hist, col = rgb(1,0,0,1/4),
     main = "CAI distribution of temperature requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency",
     cex.axis = 1) 
plot(meso.hist, col = rgb(0,0,1,1/4), add = TRUE)

legend("topleft", c("Thermophile", "Mesophile"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), bty = "n", cex=1.5))


# getting graphs with extremophiles
plot(thermo.hist, col = rgb(1,0,0,1/4),
     main = "CAI distribution of temperature requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency") 
plot(meso.hist, col = rgb(0,0,1,1/4), add = TRUE)
plot(hypther.hist, col = rgb(1,1,1/4,1/2), add = TRUE)
plot(psychrotrop.hist, col = rgb(0,1,0,1/4), add = TRUE)

legend("topleft", c("Thermophile", "Mesophile", "Hyperthermophile", "Psychrotrophic"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(1,1,1/4,1/2), rgb(0,1,0,1/4), bty = "n", cex=1.5))


# pooling extremophiles
extreme <- c(hypther.cai, psychrotrop.cai, psychro.cai)
nonextreme <- c(meso.cai, thermo.cai)



# draw the histogram plot
xlim <- range(temp.df[,3], na.rm = TRUE)
breakpoints <- seq(xlim[1], xlim[2], length.out = 20)

# get the histograms and calculate percentages
par(mfrow = c(1,1))
extreme.hist <- hist(extreme,
                  breaks = breakpoints, plot = F)
extreme.hist$counts = extreme.hist$counts/sum(extreme.hist$counts)

nonextreme.hist <- hist(nonextreme,
                   breaks = breakpoints, plot = F)
nonextreme.hist$counts = nonextreme.hist$counts/sum(nonextreme.hist$counts)


plot(extreme.hist, col = rgb(1,0,0,1/4),
     main = "CAI distribution of temperature requirements",
     xlab = "mean CAI",
     ylab = "Bacterial frequency") 
plot(nonextreme.hist, col = rgb(0,0,1,1/4), add = TRUE)

legend("topleft", c("extremophiles", "non-extremophile"), pch = 15,
       col = c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), bty = "n", cex=1.5))

legend("topleft", c("Aerobe", "Anearobe", "Facultative", "Microareophilic", 
                    "Obligate aerobe", "Obligate anaerobe"), pch = 15,
       col = c(rgb(1,1,1/4,1/2), rgb(0,0,1,1/4), rgb(1,0,0,1/4), 
               rgb(0,1,0,1/4), rgb(0,1,1,1/2), rgb(1,1/4,1,1/2)), bty = "n", cex=1.5)

