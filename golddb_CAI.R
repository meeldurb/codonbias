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


salinity <- as.data.frame(c(gold.data[1], gold.data[48]))
salinity <- na.omit(salinity)
