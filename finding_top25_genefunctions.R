#!usr/bin/env Rscript


#############################################################################
####  Author: Melanie van den Bosch
####  Title: finding genefunctions of top 25 genes 
####  Purpose of script: retrieving the genefunctions of the top 25 genes that 
####  were obtained from the iterated weight tables script output
#################################################################################



setwd("~/Documents/Master_Thesis_SSB/git_scripts")


# get filenames
outfolder <- "top25genes/"  
if (!file.exists(outfolder))dir.create(outfolder)
genomeID <- "GCA_000014005"

genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char


for (genomeID in genome.and.organisms[,1]) {
  cat (genomeID, "\n")
  fileout <- paste(outfolder, genomeID, "_top25.csv", sep="")
  if (!file.exists(fileout)) { 
    top25.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_restop25.csv", sep = "")
    genenames.files <- paste("CDS_data_withgenenames/", genomeID, "_CDS_names.csv", sep = "")
    if (file.exists(top25.files)){
      if (file.exists(genenames.files)){
        # read the top 25 file
        top25.data <- read.csv(file = top25.files, sep = ",", header = FALSE, as.is = TRUE)
        # read the CDS with names file
        genenames.data <- read.csv(file=genenames.files, sep = ",", header = TRUE)
        # find the top 25 genes in the names file
        top25 <- as.vector(top25.data[,1])
        match <- genenames.data[genenames.data[,1] %in% top25, ]
      
        # link both and write new file
        write.table(match, file=fileout, append=F, sep = ",",
                  row.names = F, quote=F, col.names=T )
      }
    }
  }
}
# 
# # read the top 25 file
# top25.files <- paste("Iterated_weight_tables_ENA/", genomeID, "_restop25.csv", sep = "")
# top25.data <- read.csv(file = top25.files, sep = ",", header = FALSE, as.is = TRUE)
# 
# # read the CDS with names file
# genenames.files <- paste("CDS_data_withgenenames/", genomeID, "_CDS_names.csv", sep = "")
# genenames.data <- read.csv(file=genenames.files, sep = ",", colClasses=c(NA, NA, NA, "NULL"), header = TRUE)
# 
# 
# # find the top 25 genes in the names file
# top25 <- as.vector(top25.data[,1])
# match <- genenames.data[genenames.data[,1] %in% top25, ]
# 
# # link both and write new file
# write.table(match, file=fileout, append=F, sep = ",",
#             row.names = F, quote=F, col.names=T )
