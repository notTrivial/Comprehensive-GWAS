# GWAS_Pipeline
GWAS Pipeline with GAPIT & TASSEL
Created by Gabrielle Dagasso
Version 0.4
Date Created: May 5, 2020
Date Updated: June 1, 2020


############################################################
This pipeline runs through:
-Structure
-Evanno Method to determine best K Value 
-GAPIT with GLM & MLM
-TASSEL with GLM & MLM

Installation of GAPIT, TASSEL, and Structure is necessary, download the pipeline folder and store softwares in it as pipeline is set to take and receive data from previous software results in appropriately named folders. 

This pipeline enables the user to determine the correct population structure K-Value with the use of Structure & Evanno Method then uses two GWAS softwares to find significant SNP's. All SNP's are adjusted using the FDR p-value adjustment procedure, and results are merged from the GAPIT/TASSEL results and outputted to the results folder

For structure you must input the max value of K you wish to test for, and for each K structure repeats 10 times so as to have better results when using the Evanno method.


############################################################
Running the Pipeline

Takes in Arguments as follows
Arg 1: is numeric genotype file with no marker names or headings
Arg 2: is the normal genotype file
Arg 3: is the phenotypic file
Arg 4: is the Kinship matrix file
Arg 5: The tassel formatted phenotype
Arg 6: the kinship matrix in proper TASSEL format (need number of individuals on first line)


Example Command
Rscript Pipeline_May30.R ./Data/Greenhouse_Numeric_Genotype ./Data/genotype_GreenHouse_editedMar1.hmp.txt ./Data/greenHouse2020_editedMissing.txt ./Data/Kinship_GreenHouse_Mar1.txt ./Data/greenHouse2020_editedMissing_tassel.txt ./Data/Kinship_GreenHouse_Mar1_tassel.txt 

############################################################
References

-Alexander E. Lipka, Feng Tian, Qishan Wang, Jason Peiffer, Meng Li, Peter J. Bradbury, Michael A. Gore, Edward S. Buckler, Zhiwu Zhang, GAPIT: genome association and prediction integrated tool, Bioinformatics, Volume 28, Issue 18, 15 September 2012, Pages 2397–2399, https://doi.org/10.1093/bioinformatics/bts444
-Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.
-Porras-Hurtado L, Ruiz Y, Santos C, Phillips C, Carracedo A, Lareu MV. An overview of STRUCTURE: applications, parameter settings, and supporting software. Front Genet. 2013;4:98. Published 2013 May 29. doi:10.3389/fgene.2013.00098

