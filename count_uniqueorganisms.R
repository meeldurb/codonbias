genome.and.organisms <- read.csv(file = "genomes_ENA.csv", header = FALSE, 
                                 as.is=TRUE) #as.is to keep the it as char



only.orgname <- sub("(\\S*\\s+\\S+) .*", '\\1', genome.and.organisms[,2])
table(only.orgname)
length(unique(only.orgname))


NeutGen <- read.csv(file= "GCneutralgenomes.csv", header = FALSE,
                    as.is = TRUE)
genome.and.organisms2 <- genome.and.organisms[which(NeutGen[,1] %in% genome.and.organisms[,1]),1:2]

only.orgname <- sub("(\\S*\\s+\\S+) .*", '\\1', genome.and.organisms2[,2])
table(only.orgname)
length(unique(only.orgname))
orgname.occur <- rle(sort(only.orgname))
morethan10 <- orgname.occur$values[which(orgname.occur$lengths > 5)]


botzman.data <- read.table(file = "botzman_list.csv", sep="," , header = TRUE,
                           as.is=TRUE)

only.orgname <- sub("(\\S*\\s+\\S+) .*", '\\1', botzman.data[,1])
table(only.orgname)
length(unique(only.orgname))
orgname.occur <- rle(sort(only.orgname))
morethan10 <- orgname.occur$values[which(orgname.occur$lengths > 5)]


git domain.occurence <- rle(sort(data[,1]))
# retrieving CAI values and mean CAI values of duplicated and unique domains
unique.domains <- domain.occurence$values[which(domain.occurence$lengths==1)]
