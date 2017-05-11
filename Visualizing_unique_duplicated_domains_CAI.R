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
library("fBasics")
library("pwr")
library("effsize")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")

genomeID <- "GCA_000003645"
cai.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)



#####________difference between unique and duplicated domains, without sampling________#####


genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

length(genome.and.organisms[,1])

# make empty columns to later bind them to dataframe
genomeIDcol <- NULL
unique.col <- NULL
dupl.col <- NULL
pvec <- numeric()

for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
  if (file.exists(cai.files)){
    data <- read.csv(file = cai.files, sep = ",", 
                     header = FALSE, as.is = TRUE)
    # how many times domains occur
    domain.occurence <- rle(sort(data[,1]))
    # retrieving CAI values and mean CAI values of duplicated and unique domains
    unique.domains <- domain.occurence$values[which(domain.occurence$lengths==1)]
    unique.domains.cai <- data[which(data[,1]%in% unique.domains),2]
    unique.domains.mean <- mean(unique.domains.cai)
    
    duplicated.domains <- domain.occurence$values[which(domain.occurence$lengths!=1)]
    duplicated.domains.cai <- data[which(data[,1]%in% duplicated.domains),2]
    duplicated.domains.mean <- mean(duplicated.domains.cai)
    sd(duplicated.domains.cai)
    sd(unique.domains.cai)
    
    # fill the columns with each iteration
    genomeIDcol <- c(genomeIDcol, genomeID)
    unique.col <- c(unique.col, unique.domains.mean)
    dupl.col <- c(dupl.col, duplicated.domains.mean)
    
    # Do a t-test to get the pvalues
    ttest <- t.test(unique.domains.cai, duplicated.domains.cai)
    pval <- ttest$p.value
    # make a vector with all the pvalues
    pvec <- c(pvec, pval)
    # by calculating the p adjusted values we correct for multiple hypothesis testing
  }
}
padj <- round(p.adjust(pvec, method = "BH", n = length(genome.and.organisms[,1])), 4)

# bind the filled columns into a dataframe
unique.duplicated.df <- data.frame(genomeIDcol, unique.col, dupl.col,
                                   padj, stringsAsFactors = FALSE)

### add factor yes/no to the df 
unique.duplicated.df$significant <- ifelse(unique.duplicated.df$padj > 0.05, "No", "Yes")
unique.duplicated.df$significant <- as.factor(unique.duplicated.df$significant)



save(unique.duplicated.df, file = "UniqueDuplicatedPadjMeanCAI.RData")


########____________ Results and draw the graph with using the p adjusted values ___________#######
load("UniqueDuplicatedPadjMeanCAI.RData")


print(paste(length(which(unique.duplicated.df[,5] == "Yes" )), "out of", 
            length(unique.duplicated.df[,5]),
            "samples are found to have a significant difference"))


unique.duplicated.df <- unique.duplicated.df[, c(2,3,5)]



myplot <- ggplot(unique.duplicated.df, aes(unique.duplicated.df[,1], 
                                           unique.duplicated.df[,2], 
                                        color=unique.duplicated.df[,3]))+ #these commands creat the plot, but nothing appears
  geom_point(size = 2, shape = 18) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  theme(legend.text = element_text(size = 15))  +
  xlab("mean CAI unique protein domains") +   
  ylab("mean CAI duplicated protein domains") +
  #ggtitle(title) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  scale_colour_discrete(name = "Significant")
#scale_fill_brewer(palette = "Spectral")

myplot   #have a look at the plot 



#####________difference between unique and duplicated domains, with sampling the data________#####

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

# set empty columns, these are filled in each iteration of the for loop
pvalcol <- NULL
genomeIDcol <- NULL
unique.col <- NULL
duplicated.col <- NULL
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")
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
    
    # the sample sizes are too different, and because a lot of observations we will get significant 
    # results alway. We therefore sample the unique and duplicated domains to reduce the power of the test
    # and also increase the effect size. In this way, the significant difference is actually meaningfull
    samplesize=500
    sampled.uniq.dom.cai <- sample(unique.domains.cai, size = samplesize, replace = TRUE)
    sampled.dupl.dom.cai <- sample(duplicated.domains.cai, size = samplesize, replace = TRUE)
    
    # calculate means of sampled domains
    unique.sampled.mean <- mean(sampled.uniq.dom.cai)
    duplicated.sampled.mean <- mean(sampled.dupl.dom.cai)
    sd(sampled.uniq.dom.cai)
    sd(sampled.dupl.dom.cai)
    # fill the columns with each iteration
    unique.col <- c(unique.col, unique.sampled.mean)
    duplicated.col <- c(duplicated.col, duplicated.sampled.mean)
    genomeIDcol <- c(genomeIDcol, genomeID)
    
    # Do a t-test to get the pvalues
    ttest <- t.test(sampled.uniq.dom.cai, sampled.dupl.dom.cai)
    pval <- ttest$p.value
    # make a vector with all the pvalues
    pvalcol <- c(pvalcol, pval)
    # by calculating the p adjusted values we correct for multiple hypothesis testing
  }
}

padj <- round(p.adjust(pvalcol, method = "BH", n = length(genome.and.organisms[,1])), 4)


# bind the filled columns into a dataframe
sampled.unique.duplicated.df <- data.frame(genomeIDcol, unique.col, dupl.col,
                                   padj, stringsAsFactors = FALSE)

### add factor yes/no to the df 
sampled.unique.duplicated.df$significant <- ifelse(sampled.unique.duplicated.df$padj > 0.05, "No", "Yes")
sampled.unique.duplicated.df$significant <- as.factor(sampled.unique.duplicated.df$significant)



save(sampled.unique.duplicated.df, file = "SampledUniqueDuplicatedPadjSign.RData")


########____________ Results and draw the graph with sampled data ___________#######


load("SampledUniqueDuplicatedPadjSign.RData")

print(paste(length(which(sampled.unique.duplicated.df[,5] == "Yes" )), "out of", 
            length(sampled.unique.duplicated.df[,5]), 
            "samples are found to have a significant difference"))


sampled.unique.duplicated.df <- sampled.unique.duplicated.df[, c(2,3,5)]



myplot <- ggplot(sampled.unique.duplicated.df, aes(sampled.unique.duplicated.df[,1], 
                                           sampled.unique.duplicated.df[,2], 
                                           color=sampled.unique.duplicated.df[,3]))+ #these commands creat the plot, but nothing appears
  geom_point(size = 2, shape = 18) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  theme(legend.text = element_text(size = 15))  +
  xlab("mean CAI unique protein domains") +   
  ylab("mean CAI duplicated protein domains") +
  #ggtitle(title) +
  guides(color = guide_legend(override.aes = list(size=6))) +
  scale_colour_discrete(name = "Significant")
#scale_fill_brewer(palette = "Spectral")

myplot   #have a look at the plot  



###__________________________for 1 genome__________________________###

genomeID <- "GCA_000007105"
cai.files <- paste("new_CAI_intradomains_ENA/", genomeID, "_intradom_CAI.csv", sep = "")

data <- read.csv(file = cai.files, sep = ",", header = FALSE, as.is = TRUE)


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


# calculating characteristics of this data
# do a t-test and check the effect size and power of the tests
unique.mean <- mean(unique.domains.cai)
unique.mean
mean(data[!duplicated(data[,1]),2]) # looking for difference
length(unique.domains.cai)
sd(unique.domains.cai)
skewness(unique.domains.cai)
duplicated.mean <- mean(duplicated.domains.cai)
duplicated.mean
mean(data[duplicated(data[,1]),2]) # looking for difference
length(duplicated.domains.cai)
sd(duplicated.domains.cai)
skewness(duplicated.domains.cai)
t.test(unique.domains.cai, duplicated.domains.cai)
# calculating effectsize and power of the test
effsize <- cohen.d(unique.domains.cai, duplicated.domains.cai, pooled = F, paired = F, conf.level = 0.95)
effsize$estimate
pwr.t.test(samplesize, d = effsize$estimate, sig.level = 0.01,
           alternative = "two.sided")

plot(duplicated.mean, unique.mean, type = "p", add = TRUE, main = "Mean CAI values of duplicated and unique protein domains",
     xlab = "Mean CAI value duplicated domain", ylab = "Mean CAI value unique domains", col = rgb(0,0,1,1/4))



# sampling to increase effect sizes and decrease power
samplesize=500
# the sample sizes are too different, we try to sample the datasets
sampled.unique <- sample(unique.domains.cai, samplesize, replace = TRUE)
sampled.duplicated <- sample(duplicated.domains.cai, samplesize, replace = TRUE)
# calculating characteristics of this data
# do a t-test and check the effect size and power of the tests
s.unique.mean <- mean(sampled.unique)
s.unique.mean
sd(sampled.unique)
s.duplicated.mean <- mean(sampled.duplicated)
s.duplicated.mean
sd(sampled.duplicated)
t.test(sampled.unique, sampled.duplicated)
# calculating effectsize and power of the test
effsize <- cohen.d(sampled.unique, sampled.duplicated, pooled = F, paired = F, conf.level = 0.95)
effsize$estimate
pwr.t.test(samplesize, d = effsize$estimate, sig.level = 0.01,
           alternative = "two.sided")


# setting the properties for drawing the graphs
xlim <- range(data[,2], na.rm = TRUE)
nbins = 25
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
for (i in 1:n.iterations){
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
for (i in 1:n.iterations){
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
