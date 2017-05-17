#!usr/bin/env Rscript

###################################################################
## Author: Melanie van den Bosch
## Script for computing the CAI of the catalase gene 
## and observing whether it is higher than the mean CAI of a genome
###################################################################

# Packages that need to be installed
source("http://bioconductor.org/biocLite.R")
?BiocUpgrade
biocLite("Biostrings")

# loading required libraries
library("Biostrings")
library("seqinr")

source("/home/melanie/Documents/Master_Thesis_SSB/git_scripts/cCAI.R")

setwd("~/Documents/Master_Thesis_SSB/git_scripts")



genomeID <- "GCA_000005825"
# reading a .csv file containing the genome names in the first column
genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char

outfolder <- "Catalase_genes/"  
if (!file.exists(outfolder))dir.create(outfolder)


####__________________after catalase seq files were written, calculate CAI and draw plot______________________#####

genomecount = 0
genome.col <- NULL
meancai.col <- NULL
catacai.col <- NULL
for (genomeID in genome.and.organisms[,1]){
  cat (genomeID, "\n")
  cai.files <- paste("new_CAI_CDS/", genomeID, "_CAI_CDS_new.csv", sep = "")
  if (file.exists(cai.files)){
    cai.data <- read.csv(file = cai.files, header = TRUE, 
                         as.is=TRUE) #as.is to keep the it as char
    mean.cai <- mean(cai.data[,2])
    
    w.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_it_weight.csv", sep = "")
    if (file.exists(w.files)){
      w.data <- read.csv(file = w.files, header = FALSE, as.is = TRUE)
      # ordering on rownames, for it to be correctly used by cai function
      ordered.w <- w.data[with(w.data, order(V1)), ]
      # only leaving numbers
      w <- ordered.w[,2]
      
      
      filein <- paste(outfolder, genomeID, "_cat.csv", sep="")
      if (file.exists(filein)){
        
        catalase.seqs <- read.csv(file = filein, sep = ",", header = TRUE, as.is = TRUE)
        # computing the postions of the domains with their corresponding postion in the CDS
        for (row in catalase.seqs){
          CDS_begin <- 3*(catalase.seqs[5]-1)+1
          CDS_end <- 3*(catalase.seqs[6])
        }
        
        # binding these columns to already existing dataframe
        # and converting to integers for later calculation
        catalase.seqs$CDS_domain_beginpos <- as.integer(CDS_begin[,1])
        catalase.seqs$CDS_domain_end <- as.integer(CDS_end[,1])
        
        # Substracting the CDS part that is associated to the protein domain
        # from columns 5 and 6 the positions of the begin and end of the 
        # CDS associated to protein domain are derived
        # when CDS are equal, take the integer from that column
        for (CDS in catalase.seqs[4]){
          dom_seq <- substr(CDS, catalase.seqs[catalase.seqs[,4] == CDS, 7], 
                            catalase.seqs[catalase.seqs[,4] == CDS, 8])
        }
        catalase.seqs$CDS_domain <- dom_seq 
        
        
        # shrinking the data to reduce computational time
        colnames(catalase.seqs) <- c("domain_ID", "domain_description", "genename", "CDS", 
                                     "d_begin","d_end", "CDS_beginpos", "CDS_end", "CDS_domain")
        keep <- c("domain_ID", "CDS_domain")
        catalase.seqs.CDS <- catalase.seqs[keep]
        
        cai.output <- compute.cai(catalase.seqs.CDS, genomeID, w, "tmpcatadomain.csv", "tmpcatadomain.fasta")
        catalase.cai <- mean(cai.output[,2])
        # fill the columns with genomeID, mean cai and mean catalase cai to later form df
        genome.col <- c(genome.col, genomeID)
        meancai.col <- c(meancai.col, mean.cai)
        catacai.col <- c(catacai.col, catalase.cai)
        
      }
    }
  }
}

catalase.df <- data.frame(genome.col, meancai.col, catacai.col,
                          stringsAsFactors = FALSE)


save(catalase.df, file = "catalase_cai.RData")


########____________ Results and draw the graph with sampled data ___________#######


load("catalase_cai.RData")


# decrease dataframe so it can be plotted by ggplot


myplot <- ggplot(catalase.df, aes(x=catalase.df[,2], 
          y=catalase.df[,3] ))+ #these commands creat the plot, but nothing appears
  geom_point(size = 2, shape = 18, color = "dodgerblue2") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 15) +
  #theme(legend.background = element_rect(fill = "white", size = .0, linetype = "dotted")) 
  #theme(legend.text = element_text(size = 15))  +
  xlab("mean CAI all genes") +   
  ylab("catalase CAI") +
  theme(legend.position="none")
  #ggtitle(title) +
#scale_fill_brewer(palette = "Spectral")

myplot   #have a look at the plot 




xlim = c(0.1, 0.9)
ylim = c(0.1, 0.9)
col = "blue"
if (genomecount == 0){  
  plot(catalase.cai, mean.cai, col=col, 
       xlim = xlim, ylim = ylim,
       pch = 18, main = "CAI complete genome vs. CAI catalase genes",
       xlab = "catalase CAI", 
       ylab = "Mean CAI")
  grid(NULL, NULL, lty = 6, col = "cornsilk2")
  #legend("bottomright" ,c("dnaE1", "dnaE2/dnaE1", "polC/dnaE3", "Other"), cex=1.5, pch=18,
  #    col=c("red", "blue", "green", "grey") , bty="n")
  genomecount = genomecount + 1
} else {
  points(catalase.cai, mean.cai, col=col, pch = 18)
}
abline(a=0, b=1)
