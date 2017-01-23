#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Plotting CAI results 
####  Purpose of script: Commands to visualize the 
####  results obtained from computing the CAI for multiple genomes
#################################################################################

install.packages("fBasics", repos="http://cran.rstudio.com/")
library("fBasics")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# intra = inside, within
# inter = between/among groups


genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


#####________Mean CAI of inter and intra domains________#####

genomecount = 0
significant = 0
nonsignificant = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.intra.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  cai.inter.files <- paste("CAI_inter_domains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    # read data and omit NA's
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    data.intra <- na.omit(data.intra)
    data.inter <- na.omit(data.inter)
    
    
    intra.cai <- mean(data.intra[,2])
    inter.cai <- mean(data.inter[,2])
    ttest <- t.test(data.intra[,2], data.inter[,2])
    pval <- ttest$p.value
    xlim <- c(0.5, 1)
    if (pval < 0.05){
      col.plot = 'blue'
      significant = significant + 1
    } 
    else {
      col.plot = 'red'
      nonsignificant = nonsignificant + 1
    }
    if (genomecount == 0){
      plot(intra.cai, inter.cai, type = 'p', xlim = xlim, ylim = xlim,
           main = "Mean CAI values of inter and intra domains",
           xlab = "Mean CAI intra domains",
           ylab = "Mean CAI inter domains",
           col = col.plot)
      genomecount = genomecount + 1
    }
    else {
      points(intra.cai, inter.cai, col=col.plot)
      legend ('topright', c('Significant', 'Non-significant'), cex = 1.5, pch = 1,
              col = c('blue', 'red'), bty='n')
    }
  }
}

total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to have a significant difference "))


#####_______Mean CAI of inter and intra domains, with sampling the data________#####



genomecount = 0
significant = 0
nonsignificant = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.intra.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  cai.inter.files <- paste("CAI_inter_domains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    # read data and omit NA's
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    data.intra <- na.omit(data.intra)
    data.inter <- na.omit(data.inter)
    

    intra.cai <- mean(data.intra[,2])
    inter.cai <- mean(data.inter[,2])
    ttest <- t.test(data.intra[,2], data.inter[,2])
    # calculate means and do a t-test to check whether observed difference is significant
    
    samplesize=500
    # the sample sizes are too different, we try to sample the datasets
    sampled.intra.cai <- sample(data.intra[,2], size = samplesize, replace = TRUE)
    sampled.inter.cai <- sample(data.inter[,2], size = samplesize, replace = TRUE)
    sampled.intra.mean <- mean(sampled.intra.cai)
    sampled.inter.mean <- mean(sampled.inter.cai)
    ttest <- t.test(sampled.intra.cai, sampled.inter.cai)
    
    pval <- ttest$p.value
    xlim <- c(0.5, 1)
    if (pval < 0.05){
      col.plot = 'blue'
      significant = significant + 1
      } 
    else {
      col.plot = 'red'
      nonsignificant = nonsignificant + 1
      }
    if (genomecount == 0){
      plot(sampled.intra.mean, sampled.inter.mean, type = 'p', xlim = xlim, ylim = xlim,
           main = "Mean CAI values of inter and intra domains (sampled data)",
           xlab = "Mean CAI intra domains",
           ylab = "Mean CAI inter domains",
           col = col.plot)
      genomecount = genomecount + 1
      }
    else {
      points(sampled.intra.mean, sampled.inter.mean, col=col.plot)
      legend ('topright', c('Significant', 'Non-significant'), cex = 1.5, pch = 1,
              col = c('blue', 'red'), bty='n')
    }
  }
}

total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to have a significant difference "))


#####________difference between CAI of inter and intra domains, computed for all protein domains________#####

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

samplesize=500
genomes.sampled <- sample(genome.and.organisms[,1], size = samplesize, replace = FALSE)


genomecount = 0
significant = 0
nonsignificant = 0
for (genomeID in genomes.sampled){
  cat (genomeID, "\n")
  cai.intra.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  cai.inter.files <- paste("CAI_inter_domains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    # read data 
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    ttest <- t.test(data.intra[,2], data.inter[,2])
    pval <- ttest$p.value
    
    xlim <- c(0, 1)
    if (pval < 0.05){
      col.plot = 'blue'
      significant = significant + 1
    } 
    else {
      col.plot = 'red'
      nonsignificant = nonsignificant + 1
    }
    if (genomecount == 0){
      plot(data.intra[,2], data.inter[,2], type = 'p', xlim = xlim, ylim = xlim, 
           main = "CAI values of inter and intra domains",
           xlab = "CAI intra domains",
           ylab = "CAI inter domains",
           col = col.plot)
      genomecount = genomecount + 1
    }
    else {
      points(sampled.intra, sampled.inter, col=col.plot)
      legend ('topleft', c('Significant', 'Non-significant'), cex = 1.5, pch = 1,
              col = c('blue', 'red'), bty='n')
    }

  }
}

total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to have a significant difference "))


###__________________________for 1 genome__________________________###

genomeID <- "GCA_000003925"
cai.intra.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
cai.inter.files <- paste("CAI_inter_domains_ENA/", genomeID, "_CAI_inter.csv", sep = "")

data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
data.intra <- na.omit(data.intra)
data.inter <- na.omit(data.inter)

length(data.intra[,2])
length(data.inter[,2])
skewness(data.intra[,2])
skewness(data.inter[,2])
intra.cai <- mean(data.intra[,2])
inter.cai <- mean(data.inter[,2])
t.test(data.intra[,2], data.inter[,2])


# setting the properties for drawing the graphs
xlim <- range(data.intra[,2], data.inter[,2], na.rm = TRUE)
nbins = 50
par(mfrow = c(2,1))


# drawing histograms of distribution cai values of unique and duplicated domains
intra.hist <- hist(data.intra[,2], 
                    breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

inter.hist <- hist(data.inter[,2], 
                        breaks=seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

# Overlap unique and duplicated domains histogram

# setting y-axis limit to draw the combined histogram
ylim <- c(0,max(inter.hist$density, inter.hist$density))

par(mfrow = c(1,1))
plot(intra.hist, col = rgb(0,0,1,1/4), xlim = xlim,  freq=FALSE,
     xlab="CAI",cex = 1.5, main="CAI distribution inter and intra domains", 
     ylab = "No. of protein domains", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

plot(inter.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)

legend("topleft" ,c("intradomain", "interdomain"), pch = 15,
       col = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)) , bty = "n",cex=1.5)

# drawing boxplots of the cai values of unique and duplicated domains
boxplot(data.intra[,2],  xlim=c(0.5,2.5), outline=FALSE, col= rgb(0,0,1,1/4),
        ylab="CAI", main = "Boxplot of inter and intra domains", cex.axis=1.5, cex.lab=1.5)
boxplot(data.inter[,2], at=2 , add=TRUE , outline=FALSE, col=rgb(1,0,0,1/4),
        cex.axis=1.5)

legend("top" ,c("Intradomains", "Interdomains"), pch=15,
       col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)) , bty="n",cex=1.5)



