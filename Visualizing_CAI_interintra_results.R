#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: Plotting CAI results 
####  Purpose of script: Commands to visualize the 
####  results obtained from computing the CAI for multiple genomes
#################################################################################

install.packages("fBasics", repos="http://cran.rstudio.com/")
install.packages("pwr", repos="http://cran.rstudio.com/")
install.packages("effsize", repos="http://cran.rstudio.com/")
install.packages("ggplot2", repos="http://cran.rstudio.com/")
library("fBasics")
library("pwr")
library("effsize")
library(ggplot2)



setwd("~/Documents/Master_Thesis_SSB/git_scripts")

# intra = inside, within
# inter = between/among groups

# open and read file with genomeIDs
genomeID <- "GCA_000003925"
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


#####__________Make df with inter intra results________####


# make empty columns to later bind them to dataframe
genomeIDcol <- NULL
intercol <- NULL
intracol <- NULL
effcol <- NULL
pwrcol <- NULL
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    if (file.exists(cai.intra.files)){
      # read data and omit NA's
      data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
      data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
      data.intra <- na.omit(data.intra)
      data.inter <- na.omit(data.inter)
      # get mean CAI of protein domains and part in betweeen
      intra.mean <- mean(data.intra[,2])
      inter.mean <- mean(data.inter[,2])
      sd(data.intra[,2])
      sd(data.inter[,2])
      
      # calculating effectsize and power of the test
      samplesize <- length(data.intra[,1])
      effsize <- cohen.d(data.intra[,2], data.inter[,2], pooled = F, paired = F, conf.level = 0.95)
      effsize$estimate
      pwr <- pwr.t.test(samplesize, d = effsize$estimate, sig.level = 0.05,
                        alternative = "two.sided")
      
      # fill columns of eff size and power
      effcol <- c(effcol, effsize$estimate)
      pwrcol <- c(pwrcol, pwr$power)
      # fill the columns by each iteration
      genomeIDcol <- c(genomeIDcol, genomeID)
      intracol <- c(intracol, intra.mean)
      intercol <- c(intercol, inter.mean)
    }
  }
}
# bind the filled columns into a dataframe
interintradf <- data.frame(genomeIDcol, intracol, intercol,
                           stringsAsFactors = FALSE)

interintradf2 <- data.frame(genomeIDcol, intracol, intercol,
                            effcol, pwrcol, stringsAsFactors = FALSE)

save(interintradf, file = "InterIntraMeans.RData")
load("InterIntraMeans.RData")


#####________calculate p-adjusted values for inter/inter domains________#####

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char
genomeID <- "GCA_000008205"

# first run the loop to create a vector with p values
pvec = numeric()
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    if (file.exists(cai.intra.files)){
    # read data and omit NA's
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    data.intra <- na.omit(data.intra)
    data.inter <- na.omit(data.inter)
    
    
    intra.cai <- mean(data.intra[,2])
    inter.cai <- mean(data.inter[,2])
    sd(data.intra[,2])
    sd(data.inter[,2])
    # Do a t-test to get the pvalues
    ttest <- t.test(data.intra[,2], data.inter[,2])
    pval <- ttest$p.value
    # make a vector with all the pvalues
    pvec <- c(pvec, pval)
    # by calculating the p adjusted values we correct for multiple hypothesis testing
    padj <- round(p.adjust(pvec, method = "BH", n = length(genome.and.organisms[,1])), 4)
    }
  }
}
   
### bind the p adjusted values to the dataframe
padj.interintradf <- cbind(interintradf, padj)
save(padj.interintradf, file = "InterIntraMeansAndPadj.RData")

### add factor yes/no to the df 
padj.interintradf$significant <- ifelse(padj.interintradf$padj > 0.05, "No", "Yes")
padj.interintradf$significant <- as.factor(padj.interintradf$significant)
save(padj.interintradf, file = "InterIntraMeansPadjSign.RData")


########____________ Results and draw the graph with using the p adjusted values ___________#######
load("InterIntraMeansPadjSign.RData")

color <- palette(c("brown3", "green4"))
palette(c("brown3", "green4"))
xlim <- c(0.3, 1)
ylim <- c(0.3, 1)

plot(padj.interintradf[,2], padj.interintradf[,3], pch = 18, 
     xlim = xlim, ylim = ylim,
     #main = "Mean CAI values of inter and intra domains",
     xlab = "Mean CAI intra domains",
     ylab = "Mean CAI inter domains",
     cex.axis = 1,
     cex.lab = 1.5,
     col = as.factor(padj.interintradf[,5]))

legend ('bottomright', c('Non-significant', 'Significant'), cex = 2, pch = 18,
        col = color, bty='n')
abline(a=0, b=1)

print(paste(length(which(padj.interintradf[,5] == "Yes" )), "out of", 
            length(padj.interintradf[,5]), 
            "samples are found to have a significant difference"))


# decrease dataframe so it can be plotted by ggplot
padj.interintradf <- padj.interintradf[, c(2,3,5)]

myplot <- ggplot(padj.interintradf, aes(padj.interintradf[,1], 
                                        padj.interintradf[,2], 
                                        color=padj.interintradf[,3]))+ #these commands creat the plot, but nothing appears
  geom_point(size = 2, shape = 18) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  theme(legend.text = element_text(size = 15))  +
  xlab("mean CAI protein domains") +   
  ylab("mean CAI non-domains") +
  #ggtitle(title) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  scale_colour_discrete(name = "Significant")
  #scale_fill_brewer(palette = "Spectral")

myplot   #have a look at the plot 



#####_______Mean CAI of inter and intra domains, with sampling the data________#####


genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# set empty columns, these are filled in each iteration of the for loop
pvalcol <- NULL
genomeIDcol <- NULL
intercol <- NULL
intracol <- NULL
#padj <- NULL
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    if (file.exists(cai.intra.files)){
      
    # read data and omit NA's
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    data.intra <- na.omit(data.intra)
    data.inter <- na.omit(data.inter)
    
    # the sample sizes are too different, and because a lot of observations we will get significant 
    # results alway. We therefore sample the unique and duplicated domains to reduce the power of the test
    # and also increase the effect size. In this way, the significant difference is actually meaningfull
    samplesize=500
    sampled.intra.cai <- sample(data.intra[,2], size = samplesize, replace = TRUE)
    sampled.inter.cai <- sample(data.inter[,2], size = samplesize, replace = TRUE)
    sd(sampled.intra.cai)
    sd(sampled.inter.cai)
    sampled.intra.mean <- mean(sampled.intra.cai)
    sampled.inter.mean <- mean(sampled.inter.cai)
    
    # fill the columns by each iteration
    genomeIDcol <- c(genomeIDcol, genomeID)
    intracol <- c(intracol, sampled.intra.mean)
    intercol <- c(intercol, sampled.inter.mean)
    
    # Do a t-test to get the pvalues
    ttest <- t.test(sampled.intra.cai, sampled.inter.cai)
    pval <- ttest$p.value
    # make a vector with all the pvalues
    pvalcol <- c(pvalcol, pval)
    # by calculating the p adjusted values we correct for multiple hypothesis testing
    }
  }
}

padj <- round(p.adjust(pvalcol, method = "BH", n = length(genome.and.organisms[,1])), 4)


# bind the filled columns into a dataframe
sampled.interintradf <- data.frame(genomeIDcol, intracol, intercol, padj,
                           stringsAsFactors = FALSE)
### add factor yes/no to the df 
sampled.interintradf$significant <- ifelse(sampled.interintradf$padj > 0.05, "No", "Yes")
sampled.interintradf$significant <- as.factor(sampled.interintradf$significant)

save(sampled.interintradf, file = "InterIntraSampledPadj.RData")

    
########____________ Results and draw the graph with sampled data ___________#######


load("InterIntraSampledPadj.RData")

print(paste(length(which(sampled.interintradf[,5] == "Yes" )), "out of", 
            length(sampled.interintradf[,5]), 
            "samples are found to have a significant difference"))


# decrease dataframe so it can be plotted by ggplot
sampled.interintradf <- sampled.interintradf[, c(2,3,5)]



myplot <- ggplot(sampled.interintradf, aes(sampled.interintradf[,1], 
                                        sampled.interintradf[,2], 
                                        color=sampled.interintradf[,3]))+ #these commands creat the plot, but nothing appears
  geom_point(size = 2, shape = 18) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  theme(legend.text = element_text(size = 15))  +
  xlab("mean CAI protein domains") +   
  ylab("mean CAI non-domains") +
  #ggtitle(title) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  scale_colour_discrete(name = "Significant")
#scale_fill_brewer(palette = "Spectral")

myplot   #have a look at the plot 




#####________difference between CAI of inter and intra domains, computed for all protein domains________#####

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

samplesize=500
genomes.sampled <- sample(genome.and.organisms[,1], size = samplesize, replace = FALSE)


pvec = numeric()
for (genomeID in genomes.sampled){
  cat (genomeID, "\n")
  cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files) && file.exists(cai.intra.files)){
    # read data 
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
    # Do a t-test to get the pvalues
    ttest <- t.test(data.intra[,2], data.inter[,2])
    pval <- ttest$p.value
    # make a vector with all the pvalues
    pvec <- c(pvec, pval)
    # by calculating the p adjusted values we correct for multiple hypothesis testing
    padj <- round(p.adjust(pvec, method = "BH", n = length(genome.and.organisms[,1])), 4)
  }
}

# then loop again over the genomes to calculate the significance by using the adjusted p-values

genomecount = 0
significant = 0
nonsignificant = 0
pvalcount = 1
for (genomeID in genomes.sampled){
  cat (genomeID, "\n")
  cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")
  if (file.exists(cai.inter.files)){
    # read data 
    data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
    data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)

    # go through all the p adjusted values, by using an iterator it will go through the padj vector and link these
    # values to the significancy of a genome
    if (padj[pvalcount] < 0.05){
      col.plot = 'blue'
      significant = significant + 1
      pvalcount = pvalcount + 1
    } 
    else {
      col.plot = 'red'
      nonsignificant = nonsignificant + 1
      pvalcount = pvalcount + 1
      }
    xlim <- c(0, 1)
    if (genomecount == 0){
      plot(data.intra[,2], data.inter[,2], type = 'p', xlim = xlim, ylim = xlim, 
           main = "CAI values of inter and intra domains",
           xlab = "mean CAI intra domains",
           ylab = "mean CAI non-domains,
           col = col.plot)
      genomecount = genomecount + 1
    }
    else {
      points(data.intra[,2], data.inter[,2], col=col.plot)
      legend ('topleft', c('Significant', 'Non-significant'), cex = 1.5, pch = 1,
              col = c('blue', 'red'), bty='n')
    }

  }
}

total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to have a significant difference "))


###__________________________for 1 genome__________________________###

genomeID <- "GCA_000003645"
cai.intra.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
cai.inter.files <- paste("new_CAI_interdomains_ENA/", genomeID, "_CAI_inter.csv", sep = "")

data.intra <- read.csv(file = cai.intra.files, sep = ",", header = FALSE, as.is = TRUE)
data.inter <- read.csv(file = cai.inter.files, sep = ",", header = FALSE, as.is = TRUE)
data.intra <- na.omit(data.intra)
data.inter <- na.omit(data.inter)

length(data.intra[,2])
length(data.inter[,2])
skewness(data.intra[,2])
skewness(data.inter[,2])
mean(data.intra[,2])
mean(data.inter[,2])
sd(data.intra[,2])
sd(data.inter[,2])

ttest <- t.test(data.intra[,2], data.inter[,2])
# calculate means and do a t-test to check whether observed difference is significant

samplesize=500
# the sample sizes are too different, we try to sample the datasets
sampled.intra.cai <- sample(data.intra[,2], size = samplesize, replace = TRUE)
sampled.inter.cai <- sample(data.inter[,2], size = samplesize, replace = TRUE)
mean(sampled.intra.cai)
mean(sampled.inter.cai)
sd(sampled.intra.cai)
sd(sampled.inter.cai)
ttest <- t.test(sampled.intra.cai, sampled.inter.cai)
pwr.t.test(samplesize, d = 0.152, sig.level = 0.01,
                    alternative = "two.sided")
# setting the properties for drawing the graphs
xlim <- range(data.intra[,2], data.inter[,2], na.rm = TRUE)
nbins = 30
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


