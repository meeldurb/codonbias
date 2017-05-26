# Codon Usage Bias
Scripts written during my Master thesis @ systems and synthetic biology WUR

## Explanation of scripts


### Computation of weight tables & results

compute_refweighttables_ENA.R: algorithm to compute the weighttables from ENA db, seed is ribosomal protein encoding genes.

cweight.R: Script with function to calculate weight tables.

iterated_refweighttables_ENA.R: weight tables computed with an iterative algorithm, uses ribosomal protein encoding genes as initial seed

flowchart_italg.png: flowchart of the iterative algorithm.

count_itweight.py: from the tmp files produced from "iterated_refweighttables_ENA.R" parsing the iteration counts into a csv file

itcount_final.csv: all the iteration counts of all the genomes obtained from ""iterated_refweighttables_ENA.R" its tmpfile and calculated with "count_itweight.py"

robustnesscheck_itweighttables_ENA.R: Uses the same algorithm as in "iterated_refweighttables_ENA.R", but uses a random seed to compute intial weight table. Also draws calculates difference between ribosomal and random seed and draws histograms from the data.


### Weight tables results
itcount_rand_final1.csv: all the iteration counts of all the genomes obtained random seed from "robustnesscheck_itweighttables_ENA.R" its tmpfile and calculated with "count_itweight.py"

Visualizing_iterativeweight_results.R: Draws graph on how many iterations each genome went through gotten from "
itcount_final.csv".

itcountgraph.png: graph of iteration counts ribosomal protein encoding genes seed to "iterated_refweighttables_ENA.R"

it1_4genomes.csv: genomes with 1-4 iterations before iterative weight tables are written.

it5_21genomes.csv: genomes with 5-21 iterations before iterative weight tables are written.

robust_vs_ribo_seed.R: Visualizing the differences between ribosomal seed and random seed.

diffriboranditcounts.png: difference between the iteration counts of ribosomal seed ans random seed.

finding_top25_genefunctions.R: Finds the top 25 obtained from "iterated_refweighttables_ENA.R" and links the genefunctions obtained from "ENADB_get_CDS_genomes_andnames.R" to it.

write_relativeadaptiveness.R: Makes a dataframe of all genomes and all codons with their corresponding weights. 

GenomeDataSet.RData: A dataframe of all genomes and all codons with their corresponding weights. 

boxplots_relativeadapt.R: Drawing boxplots of codosn of all genomes. Takes "GenomeDataset.RData" as input

boxplotsreladapt.png: The resulting plot of the boxplots of codons.

PCA_relativeadaptiveness.R: Script calculating and drawing PCA plots from the weight tables of all the genomes. Takes "GenomeDataset.RData" as input.

thirdpos_PCA_biplotresults.R: Script for visualizing relation between loadings of the PCA, the PC's and nucleotide position on 3rd position of codon. Takes "GenomeDataset.RData" as input .

codons_lefttop_Tbias.png: from the biplot we see that the genomes that are in the lefttop of the PCA have a bias towards a T at 3rd nucleotide. These genomes are substracted and boxplots are drawn.

codons_lefttop_Cbias.png: from the biplot we see that the genomes that are in the righttop of the PCA have a bias towards a C at 3rd nucleotide. These genomes are substracted and boxplots are drawn.

codons_thirdpos.png: plots together of 3 position.

PCA_relativeadaptiveness.R:  makes this in a huge dataframe "codonGenomeDataSet.RData" gotten from "iterated_refweighttables_ENA.R" and draws PCA plot of genomes and its iterative weight table grouped by certain conditions in the "gold_gca.tsv" file.


### catalase CAI
catalaseCAI.R: Calculates, draws dataframe and draws graph of the average CAI of complete genome and the CAI of catalase domains in the genome.

catalase_cai.RData: Dataframe containing mean CAI of all genes and CAI of catalase of all genomes.

catalase_cai.png: graph of mean CAI all genes and CAI of catalase.


### Computation of codon pairs & results
cCAIpairs.R: script with funciton to calculate CAI of codonpairs.

ccodpairsweight.R: Script with function to calculate codonpair weight tables.

codonpairs_itweighttable.R: Script to compute the iterative weight tables of the codonpairs. Written to folder "codonpairs_itweight/".

write_codonpairsweight.R: writes a dataframe containing all genomes with all weights of the codon pairs. To the file "codonpairsGenomeDataSet.Rdata".

codonpairsGenomeDataset.RData: dataframe containing all genomes with all weights of the codon pairs.

PCA_codonpairs.R: Drawing PCA plots and biplots of all genomes, for different taxonomic ranks. Takes input from "codonpairsGenomeDataSet.Rdata".

boxplots_codonpairs.R: Drawing boxplots of codon pairs of all genomes, each boxplot contains all codonpairs for a unique amino acid pair. Takes input from "codonpairsGenomeDataSet.Rdata".


### GC content, CAI and POLIII
GC_content_vs_CAI: writes a dataframe of all genomes with the GC content and mean CAI of all genomes.

GCcontMeanCAI.RData: dataframe of all genomes with the GC content and mean CAI of all genomes.

CAI_vs_GCcont.png: Image CAI vs. GC content.

GC_content_vs_CAI.R: Draws a graph between correlation of meanGCcont (calc GC content from sequences in "ENADB_get_CDS_genomes.R") and meanCAI obtained from "calculate_CAI_CDS.R".

GC_content_vs_CAI_POLIII.R: writes a file of genomeID, GC content (calc GC content from sequences in "ENADB_get_CDS_genomes.R") meanCAI obtained from "calculate_CAI_CDS.R" and polIII isomer (gotten from "CAI_GCcont_POLIII_allgenomes.csv"), then it draws a graph with regression lines.

CAI_GCcont_POLIII_allgenomes.csv: list from Zhao et al. containing the genomes that they tested with their genera names and polII isomerase.


### environmental lifestyle 
CAI_vs_GCcont.png: Image CAI vs. GC content

meancaidist.png: distribution of genomes and their mean CAI value for all genes. 

count_uniqueorganisms: counting the unique organims we have in the data to justify how many organims were tested in these results.

GCextremegenomes.csv: genomes with gc content below 35% and higher than 65%.

GCneutralgenomes.csv: genomes with gc content between 35% and 65%.

gold_gca.tsv: gold database containing more information associated to the genomes in ENAdb.

golddb_CAI.R: compares environmental conditions by making use of "gold_gca.tsv".

halotol.cai.png: halotolerant and halophile genomes histogram with mean CAI for all genes.

pronooxygen.cai.png: genomes that need or do not need oxygen to survive and their mean CAI value for all genes.

pronooxygen.cai.png: mesophile and thermophile genomes and their mean CAI value for all genes.


### General 
test_genomes_ENA10.csv: test file with 10 genomes.

get_genomes_ENA.R: retrieves all the genomes + organism from ENA db.

genomes_ENA.csv: all the genomes with organism obtained from the ENA db, used in almost all R scripts.

write2fasta.py: writes the csv file gotten from the scripts where CAI is computed to a fasta file for further computation of CAI.

Write_csvtofasta: 1st file that was written to convert the csv file of domains to fasta for further computation of CAI.

Write_genecsvtofasta: converts the gene sequences from csv format to fasta format for further computation of CAI.

cCAI.R: script with function to calculate CAI.

calculate_CAI_CDS.R: Calculates the CAI of all genes with iterative weight tables of all genomes in ENA db and writes files to "new_CAI_CDS/".

ENADB_get_CDS_genomes.R: retrieves all the genomes with its complete CDS from ENA db written to CDS_data/ folder.

ENADB_get_CDS_genomes_andnames.R: retrieves all the genomes with its complete CDS and gene names from ENA db.


### Domains vs. non-domains
calculate_CAI_complete_ENA.R: Calculates CAI of CDS associated to domains.

calculate_CAI_domains_ENA.R: Calulates CAI of all domains associated to genomes using old weight tables gotten from "compute_refweighttables_ENA.R".

calculate_CAI_interdomains_neww.R: Calculates CAI of all parts in between domains using iterative weight tables.

calculate_CAI_intradomains_neww.R: Calculates CAI of all parts inside domains using iterative weight tables.

ENADB_get_domains_genomes.R: retrieves all the genomes with its associated domains and DNA sequence fom ENA db and writes domain data to folder "Domain_data_ENA/".

Visualizing_CAI_interintra_results.R: Draws graph of the results obtained from "calculate_CAI_interdomains_neww.R" and "calculate_CAI_intradomains_neww.R".

InterIntraMeans.RData: dataframe of mean CAI domain and mean CAI non-domains of all genomes.

InterIntraMeansPadjSign.RData: dataframe of mean CAI domain, mean CAI non-domains, p adjusted values (correcting for multiple comparisons) and its significance yes/no of all genomes.

InterIntraSampledPadj.RData: dataframe of mean CAI domain, mean CAI non-domains, p adjusted values (correcting for multiple comparisons) and its significance yes/no of all genomes, but now sampled for only 500 observations for domains and non-domains.

interintera1: graph of domains vs non-domains with corrected for multiple testing.


### Domains unique & duplicated
Visualizing_CAI_results.R: Draws graphs of comparison between unique and duplicated domains obtained from "calculate_CAI_domains_ENA.R".

Visualizing_unique_duplicated_domains_CAI: Draws graph of comparison between unique and duplicated domains obtained from "calculate_CAI_intradomains_neww.R".

UniqueDuplicatedPadjMeanCAI.RData: dataframe of mean CAI unique domains, mean CAI duplicated domains, p adjusted values (correcting for multiple comparisons) and its significance yes/no of all genomes.

uniqueduplicated1.png: Graph of mean CAI of unique and duplicated domains.

SampledUniqueDuplicatedPadjSign.RData: dataframe of mean CAI unique domains, mean CAI duplicated domains, p adjusted values (correcting for multiple comparisons) and its significance yes/no of all genomes, but now sampled for only 500 observations for domains and non-domains.

sampleduniqueduplicated.png: Graph of mean CAI of unique and duplicated domains, but now sampled for only 500 observations for unique domains and duplicated domains.



### Other files


EC_Pfam_calculated_associations.csv: file containing association between EC numbers and Pfam ID's

















