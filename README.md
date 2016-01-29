# ratcon

ncomm_blais_check_tasks.m

This MATLAB script imports iRno and iHsa into the RAVEN toolbox and performs the checkTasks function for 327 metabolic tasks described in Supplementary Table 4. Most tasks capture metabolic functions that are common to both rats and humans. Some tasks are meant to fail in both organisms while others are meant to fail in only one organism, consistent with known species-specific differences between rat and human metabolism.




ncomm_helper.R
R source file that defines several helper functions used in the R scripts below.

ncomm_blais_gpr_conversion.R
This script was used to generate Fig. 2a and Supplementary Fig. 3. For Supplementary Fig. 3, this R script inputs orthology annotations from multiple databases by reading in data from Supplementary Table 2 to generate rat gene-protein-reaction (GPR) relationship rules based on human GPR rules from HMR2 (after replacing Ensembl Gene Identifiers with Entrez gene identifiers). For Fig. 2a, this R script compares the complexity of gene-protein-reaction (GPR) directly between rat and human metabolic models by reading in data from Supplementary Table 3. This script inputs gene-protein-reaction (GPR) relationship information from a superset of reactions included in iRno and iHsa to compare the global distribution of GPR sizes between the rat and human models. The goal of this script was to demonstrate how GPR sizes varied between species when filtering orthology annotations for the automated conversion of iHsa to iRno as well as after manual curation. 
