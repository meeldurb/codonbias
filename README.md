# Scripts for the Master Thesis on codon usage bias

One Paragraph of project description goes here

## Explanation of scripts

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.


calculate_CAI_CDS.R: Calculates the CAI of complete CDS of all genomes in ENA db.

calculate_CAI_complete_ENA.R: Calculates CAI of CDS associated to domains

calculate_CAI+domains_ENA.R: Calulates CAI of all domains associated to genomes

calculate_CAI_interdomains_neww.R: Calculates CAI of all parts in between domains using iterative weight tables

calculate_CAI_intradomains_neww.R: Calculates CAI of all parts inside domains using iterative weight tables

cCAI.R: script with function to calculate CAI

codonGenomeDataSet.RData: R dataset containing all the weight tables of all the genomes in one frame, suited for PCA

compute_refweighttables_ENA.R: algorithm to compute the weighttables from ENA db, seed is ribosomal protein encoding genes

count_itweight.py: from the tmp files produced from "iterated_refweighttables_ENA.R" parsing the iteration counts into a csv file

cweight.R: Script with function to calculate weight tables

EC_Pfam_calculated_associations.csv: file containing association between EC numbers and Pfam ID's


# 



