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

#genomeID <- "GCA_000003645"
#cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
#data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

#####________difference between unique and duplicated domains, without sampling________#####

genomecount = 0
significant = 0
nonsignificant = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  if (file.exists(cai.files)){
    data <- read.csv(file = cai.files, sep = ",", 
                     header = FALSE, as.is = TRUE)
    # how many times domains occur
    domain.occurence <- rle(sort(data[,1]))
    # retrieving CAI values of duplicated and unique domains
    unique.domains <- domain.occurence$values[which(domain.occurence$lengths==1)]
    unique.domains.cai <- data[which(data[,1]%in% unique.domains),2]
    
    duplicated.domains <- domain.occurence$values[which(domain.occurence$lengths!=1)]
    duplicated.domains.cai <- data[which(data[,1]%in% duplicated.domains),2]
    
    # calculate means and do a t-test to check whether observed difference is significant
    unique.mean <- mean(unique.domains.cai)
    duplicated.mean <- mean(duplicated.domains.cai)
    ttest <- t.test(unique.domains.cai, duplicated.domains.cai)
    pval <- ttest$p.value
    if (pval < 0.05){
      col.plot= "blue"
      significant = significant + 1
    } else{
      col.plot = "red"
      nonsignificant = nonsignificant + 1
    }
    if (genomecount == 0){
    plot(duplicated.mean, unique.mean, type = "p", xlim = c(0.5, 1), ylim = c(0.5, 1),
         main = "Mean CAI values of duplicated and unique protein domains",
         xlab = "Mean CAI duplicated", 
         ylab = "Mean CAI unique", 
         col = col.plot)
      genomecount = genomecount + 1
    } 
    else {
      points(duplicated.mean, unique.mean, col = col.plot)
      legend("topright" ,c("Unique", "Duplicated"), cex=1.5, pch=1,
             col=c("blue", "red") , bty="n")
    } 
  }
  
}

total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to be significant"))


#####________difference between unique and duplicated domains, with sampling the data________#####
genomecount = 0
significant = 0
nonsignificant = 0
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("CAI_domains_ENA/", genomeID, "_CAI.csv", sep = "")
  if (file.exists(cai.files)){
    data <- read.csv(file = cai.files, sep = ",", 
                     header = FALSE, as.is = TRUE)
    # how many times domains occur
    domain.occurence <- rle(sort(data[,1]))
    # retrieving CAI values of duplicated and unique domains
    unique.domains <- domain.occurence$values[which(domain.occurence$lengths==1)]
    unique.domains.cai <- data[which(data[,1]%in% unique.domains),2]
    
    duplicated.domains <- domain.occurence$values[which(domain.occurence$lengths!=1)]
    duplicated.domains.cai <- data[which(data[,1]%in% duplicated.domains),2]
    
    # calculate means and do a t-test to check whether observed difference is significant
    # not sure whether I still have to do this:
    samplesize=500
    # the sample sizes are too different, we try to sample the datasets
    sampled.uniq.dom.cai <- sample(unique.domains.cai, size = samplesize, replace = TRUE)
    sampled.dupl.dom.cai <- sample(duplicated.domains.cai, size = samplesize, replace = TRUE)
    sampled.unique.mean <- mean(sampled.uniq.dom.cai)
    sampled.duplicated.mean <- mean(sampled.dupl.dom.cai)
    ttest <- t.test(sampled.uniq.dom.cai, sampled.dupl.dom.cai)
    pval <- ttest$p.value
    if (pval < 0.05){
      col.plot= "blue"
      significant = significant + 1
    } else{
      col.plot = "red"
      nonsignificant = nonsignificant + 1
    }
    if (genomecount == 0){
      plot(sampled.duplicated.mean, sampled.unique.mean, type = "p", xlim = c(0.5, 1), ylim = c(0.5, 1),
           main = "Mean CAI values of duplicated and unique protein domains (sampled data)",
           xlab = "Mean CAI duplicated", 
           ylab = "Mean CAI unique", 
           col = col.plot)
      genomecount = genomecount + 1
    } 
    else {
      points(sampled.duplicated.mean, sampled.unique.mean, col = col.plot)
      legend("topright" ,c("Unique", "Duplicated"), cex=1.5, pch=1,
             col=c("blue", "red") , bty="n")
    } 
  }

}
total = significant + nonsignificant
print(paste(significant, "out of", total, "samples are found to be significant in the sampled dataset"))



###__________________________for 1 genome__________________________###
# how many times a domain occurs
domains.occurence <- rle(sort(data[,1]))

# taking the values of the domains that are unique and duplicated
unique.domains <- domains.occurence$values[which(domains.occurence$lengths==1)]
unique.domains.cai <- data[which(data[,1]%in% unique.domains),2]

duplicated.domains <- domains.occurence$values[which(domains.occurence$lengths!=1)]
duplicated.domains.cai <- data[which(data[,1]%in% duplicated.domains),2]

# same way of calculating
#domains.occurence2 <- table(data[,1])
#threshold <- 2
#unique.domains2 <- names(domains.occurence2[domains.occurence2<threshold])
#duplicated.domains2 <- names(domains.occurence2[domains.occurence2>=threshold])


# calculate means, size sd and skewness of the samples
# do a t-test to check whether observed difference is significant
unique.mean <- mean(unique.domains.cai)
mean(data[!duplicated(data[,1]),2]) # looking for difference
length(unique.domains.cai)
sd(unique.domains.cai)
skewness(unique.domains.cai)
duplicated.mean <- mean(duplicated.domains.cai)
mean(data[duplicated(data[,1]),2]) # looking for difference
length(duplicated.domains.cai)
sd(duplicated.domains.cai)
skewness(duplicated.domains.cai)
t.test(unique.domains.cai, duplicated.domains.cai)

plot(duplicated.mean, unique.mean, type = "p", add = TRUE, main = "Mean CAI values of duplicated and unique protein domains",
     xlab = "Mean CAI value duplicated domain", ylab = "Mean CAI value unique domains", col = rgb(0,0,1,1/4))



# not sure whether I still have to do this:
samplesize=500
# the sample sizes are too different, we try to sample the datasets
sampled.uniq.dom.cai <- sample(unique.domains.cai, 500, replace = TRUE)
sampled.dupl.dom.cai <- sample(duplicated.domains.cai, 500, replace = TRUE)
sampled.unique.mean <- mean(sampled.uniq.dom.cai)
sampled.duplicated.mean <- mean(sampled.dupl.dom.cai)
t.test(sampled.uniq.dom.cai, sampled.dupl.dom.cai)


# setting the properties for drawing the graphs
xlim <- range(data[,2], na.rm = TRUE)
nbins = 50
par(mfrow = c(3,1))

# drawing histogram with all data inside 
hist(data[,2], breaks=seq(xlim[1], xlim[2], length = nbins), 
     xlim = xlim)

# drawing histograms of distribution cai values of unique and duplicated domains
unique.hist <- hist(unique.domains.cai, 
                    breaks = seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

duplicated.hist <- hist(duplicated.domains.cai, 
                        breaks=seq(xlim[1], xlim[2], length = nbins), xlim = xlim)

# Overlap unique and duplicated domains histogram

# setting y-axis limit to draw the combined histogram
ylim <- c(0,max(unique.hist$density, duplicated.hist$density))

##pdf("CAI_single_or_multiple_copies.pdf", width=14)
par(mfrow = c(1,1))
plot(unique.hist, col = rgb(0,0,1,1/4), xlim = xlim,  freq=FALSE,
     xlab="CAI",cex = 1.5, main="CAI distribution of unique and non-unique protein domains", 
     ylab = "No. of protein domains", cex.axis = 1.5, cex.lab=1.5, ylim = ylim) 

plot(duplicated.hist, col = rgb(1,0,0,1/4), xlim =  xlim, add = TRUE,  freq = FALSE)

legend("topright" ,c("Unique", "Duplicated"), pch = 15,
       col = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)) , bty = "n",cex=1.5)

# drawing boxplots of the cai values of unique and duplicated domains
boxplot(unique.domains.cai,  xlim=c(0.5,2.5), outline=FALSE, col= rgb(0,0,1,1/4),
        ylab="CAI", main = "Boxplot of unique and non-unique protein domains", cex.axis=1.5, cex.lab=1.5)
boxplot(duplicated.domains.cai, at=2 , add=TRUE , outline=FALSE, col=rgb(1,0,0,1/4),
        cex.axis=1.5)

legend("top" ,c("Unique", "Non-Unique"), pch=15,
       col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)) , bty="n",cex=1.5)
#dev.off()





####*******************Testing whether the CAI is conserved in domains or in proteins******************####
# Inter vs intra variability

# Making a list of all domainsIDs that are inside the dataset, not duplicate them
domain.list=unique(data[,1]) # do not take unique though, just all domains

# sample 2 domains from the list 
# Check to which positions the sampled Pfam IDs are found in the dataset
# sample one from these and return the CAI values
# then substract 

n.iterations <- 1000
interDist = NULL
for (i in 1:Niter){
  sel=sample(domain.list,2 )
  interDist <- c(interDist, abs(data[sample(which(data[,1]==sel[1]),1),2]-
                                  data[sample(which(data[,1]==sel[2]),1),2]))
}

# Do the same for domains that occur more than once
# sample one Pfam ID from the dataset of domains that are non-unique
# then sample 2 Pfam IDs from the original dataset and see whether they are equal
# to the sample of the non-unique domains
# then substract the values of the 
intraDist=NULL
for (i in 1:Niter){
  sel=sample(moreOnce,1 )
  sel=sample(which(data[,1]==sel),2) # returns the indices of the Pfam ID that is equal to the one in the dataset
  intraDist <- c(intraDist, abs(data[sel[1],2]- data[sel[2],2])) 
#returns the numbers inside data from the indices sel[1/2],1/2 and substracts these
# this is done 10000 times
}

#Visualizing the data
pdf("CAI_differences_between_domains.pdf")
boxplot(cbind(interDist, intraDist), outline=FALSE, names=c("Different
domains", "Same domains"), ylab = "Difference in CAI", col=c("aquamarine", "coral"))
dev.off()
t.test(interDist, intraDist)
