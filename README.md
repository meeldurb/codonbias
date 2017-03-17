# Codon Usage Bias
Scripts written during my Master thesis @ systems and synthetic biology WUR

## Explanation of scripts



calculate_CAI_CDS.R: Calculates the CAI of complete CDS of all genomes in ENA db.

calculate_CAI_complete_ENA.R: Calculates CAI of CDS associated to domains

calculate_CAI_domains_ENA.R: Calulates CAI of all domains associated to genomes using old weight tables gotten from "compute_refweighttables_ENA.R"

calculate_CAI_interdomains_neww.R: Calculates CAI of all parts in between domains using iterative weight tables

calculate_CAI_intradomains_neww.R: Calculates CAI of all parts inside domains using iterative weight tables

cCAI.R: script with function to calculate CAI

compute_refweighttables_ENA.R: algorithm to compute the weighttables from ENA db, seed is ribosomal protein encoding genes

count_itweight.py: from the tmp files produced from "iterated_refweighttables_ENA.R" parsing the iteration counts into a csv file

cweight.R: Script with function to calculate weight tables

ENADB_get_CDS_genomes.R: retrieves all the genomes with its complete CDS from ENA db

ENADB_get_CDS_genomes_andnames.R: retrieves all the genomes with its complete CDS and gene names from ENA db

ENADB_get_domains_genomes.R: retrieves all the genomes with its associated domains and DNA sequence fom ENA db

finding_top25_genefunctions.R: Finds the top 25 obtained from "iterated_refweighttables_ENA.R" and links the genefunctions obtained from "ENADB_get_CDS_genomes_andnames.R" to it.

GC_content_vs_CAI.R: Draws a graph between correlation of meanGCcont (calc GC content from sequences in "ENADB_get_CDS_genomes.R") and meanCAI obtained from "calculate_CAI_CDS.R"

get_genomes_ENA.R: retrieves all the genomes + organism from ENA db.

golddb_CAI.R: compares environmental conditions by making use of "gold_gca.tsv".

iterated_refweighttables_ENA.R: weight tables computed with an iterative algorithm, uses ribosomal protein encoding genes as initial seed

PCA_relativeadaptiveness.R:  makes this in a huge dataframe "codonGenomeDataSet.RData" gotten from "iterated_refweighttables_ENA.R" and draws PCA plot of genomes and its iterative weight table grouped by certain conditions in the "gold_gca.tsv" file

robustnesscheck_itweighttables_ENA.R: Uses the same algorithm as in "iterated_refweighttables_ENA.R", but uses a random seed to compute intial weight table

Visualizing_CAI_interintra_results.R: Draws graph of the results obtained from "calculate_CAI_interdomains_neww.R" and "calculate_CAI_intradomains_neww.R"

Visualizing_CAI_results.R: Draws graphs of comparison between unique and duplicated domains obtained from "calculate_CAI_domains_ENA.R"

Visualizing_iterativeweight_results.R: Draws graph on how many iterations each genome went through gotten from "
itcount_final.csv"

Visualizing_unique_duplicated_domains_CAI: Draws graph of comparison between unique and duplicated domains obtained from "calculate_CAI_intradomains_neww.R"

write2fasta.py: writes the csv file gotten from the scripts where CAI is computed to a fasta file for further computation of CAI

Write_csvtofasta: 1st file that was written to convert the csv file of domains to fasta for further computation of CAI

Write_genecsvtofasta: converts the gene sequences from csv format to fasta format for further computation of CAI

## Explanation of files

codonGenomeDataSet.RData: R dataset containing all the weight tables of all the genomes in one frame, suited for PCA

EC_Pfam_calculated_associations.csv: file containing association between EC numbers and Pfam ID's

genomes_ENA.csv: all the genomes with organism obtained from the ENA db

genomes_ENA_nodupl.csv: all the genomes with organism, but without duplicated genomeIDs from ENA db

gold_gca.tsv: gold database containing more information associated to the genomes in ENAdb

itcount_final.csv: all the iteration counts of all the genomes obtained from ""iterated_refweighttables_ENA.R" its tmpfile and calculated with "count_itweight.py"

test_genomes_ENA10.csv: test file with 10 genomes




