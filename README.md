# DNADamage_Phylointeracome

R scripts and resources  used for the STRINGdb analysis included in the DNA damage phylointeractome project. Data from this project will be soon published under the following title:

**DNA damage repair proteins across the tree of life publication**

This repository **is oriented to readers of the publication**. The code to generate the analysis and figures related to the functional analysis (KEGG and GO terms) and the STRINGdb networks in the paper is provided. 

The repository contains one folder and two scripts.

## 00_Rawdata

This folder contains the requiered files for the analysis:

- ProteinsPerSpecies.csv. It contains the enriched proteins for each DNA lesion.

- ...protein.links.detailed.v11.5.txt. For each species, the downloaded STRINGdb interactions (version 11.5). To get more details and download the actual files go to STRING (https://string-db.org/cgi/download?sessionId=bGJTpaVuOTAX)

## STRINGdb_Analysis_PerSpecies.R

R script containing the analysis to generate the interactions networks per species as shown in **Figure 6A** and **Supplementary Figure 8**.

## STRINGdb_Analysis_PerSpeciesLesion.R

R script containing the analysis to generate, per species and lesion, the KEGG and GO analysis, as shown in **Figure 2** and **Supplemental Figure 2** and the interactions networks.