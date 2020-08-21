GWAS_Pipeline
Created by Gabrielle Dagasso & Dr. Lingling Jin
Version 0.4
Date Created: May 5, 2020
Date Updated: June 25, 2020

######################## INSTALLATION & SETUP ###################################

Download the pipeline folder available here, the folder contains subfolders for which the pipeline takes in and outputs results. 

Requires R v3.6.1, all R packages are installed automatically when running the RScript pipeline. 

Installation of TASSEL, and Structure is necessary. Install both of these softwares within the downloaded pipeline folder as pipeline is set to run from this folder. The initial working directory is the downloaded pipeline folder, the RScript from there handles the current working directories. 

Structure - https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html
TASSEL - https://www.maizegenetics.net/tassel

For Structure you must update the mainparams file according to the numerical genotypic input data you have. To transform genotype data into numerical form one can use TASSEL with the following command

./run_pipeline.pl -h [input data] -NumericalGenotypePlugin -endPlugin -export output -exportType ReferenceProbability


########################## INFORMATION ##################################
This pipeline runs through:
-Structure
-Evanno Method to determine best K Value 
-GAPIT with GLM & MLM
-TASSEL with GLM & MLM
- GLMNET with LASSO method

This pipeline enables the user to determine the correct population structure K-Value with the use of Structure & Evanno Method then uses both GAPIT & TASSEL to find significant SNP's. All SNP's are adjusted using the FDR p-value adjustment procedure, and results are merged from the GAPIT/TASSEL results and outputted to the results folder.

For structure you must input the max value of K you wish to test for, and for each K structure repeats 10 times so as to have better results when using the Evanno method.


############################ RUNNING THE PIPELINE ################################

Takes in Arguments as follows
Arg 1: numeric genotype file with no marker names or headings (as per Structure requirements)
Arg 2: non-numeric genotype file
Arg 3: phenotypic file
Arg 4: Kinship matrix file
Arg 5: TASSEL formatted phenotype file
Arg 6: the kinship matrix in proper tassel format, need number of individuals on first line
Arg 7: Numeric Genotype with SNP/Marker names


Example Command
Rscript Pipeline_June25.R ./Data/Numeric_Geno_file ./Data/geno_file.hmp.txt ./Data/PHENO_FILE.txt ./Data/Kinship.txt ./Data/PHENO_FILE_TASSEL.txt ./Data/Kinship_Matrix.txt ./Data/numeric_genotype_WithSNPNames.txt


############################################################
References

-Alexander E. Lipka, Feng Tian, Qishan Wang, Jason Peiffer, Meng Li, Peter J. Bradbury, Michael A. Gore, Edward S. Buckler, Zhiwu Zhang, GAPIT: genome association and prediction integrated tool, Bioinformatics, Volume 28, Issue 18, 15 September 2012, Pages 2397–2399, https://doi.org/10.1093/bioinformatics/bts444
-Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES. (2007) TASSEL: Software for association mapping of complex traits in diverse samples. Bioinformatics 23:2633-2635.
-Porras-Hurtado L, Ruiz Y, Santos C, Phillips C, Carracedo A, Lareu MV. An overview of STRUCTURE: applications, parameter settings, and supporting software. Front Genet. 2013;4:98. Published 2013 May 29. doi:10.3389/fgene.2013.00098

