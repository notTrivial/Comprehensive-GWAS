#GWAS Pipeline with GAPIT
#Created by Gabrielle Dagasso
#Version 0.6
#Date Created: May 5, 2020
#Date Updated: July 1, 2020
#Update 1: Added in TASSEL
#Update 2: Added in GLMNET & Fixed Merged Results Function
#Runs through Structure, evanno method to determine the correct K, GAPIT & TASSEL with GLM and MLM
#Merges TASSEL & GAPIT results and outputs manhattan plots & stats by phenotype in Merged_Results Folder
#Runs Lasso based upon most significant SNP's found in the merged results from GAPIT & TASSEL (only if sig SNPs are found)

#Takes in Arguments as follows
#Arg 1: is numeric genotype file with no marker names or headings
#./Data/Greenhouse/Greenhouse_Numeric_Genotype
#Arg 2: is the normal genotype file
#./Data/genotype_GreenHouse_editedMar1.hmp.txt
#Arg 3: is the phenotypic file
#./Data/greenHouse2020_editedMissing.txt
#Arg 4: is the Kinship matrix file
#./Data/Kinship/Kinship_GreenHouse_Mar1.txt 
#Arg 5: The tassel formatted phenotype
#./Data/greenHouse2020_editedMissing_tassel.txt
#Arg 6: the kinship matrix in proper tassel format, need number of individuals on first line
#./Data/Kinship_GreenHouse_Mar1_tassel.txt
#Arg 7: Numeric Genotype with SNP/Marker names
#./Data/numeric_genotype_WithSNPNames.txt


#Example Command
#Rscript Pipeline_June25.R ./Data/Greenhouse_Numeric_Genotype ./Data/genotype_GreenHouse_editedMar1.hmp.txt ./Data/greenHouse2020_editedMissing.txt ./Data/Kinship_GreenHouse_Mar1.txt ./Data/greenHouse2020_editedMissing_tassel.txt ./Data/Kinship_GreenHouse_Mar1_tassel.txt ./Data/numeric_genotype_WithSNPNames.txt

###################################### FUNCTIONS UTILIZED ######################################
## helper function to check if specific packagens have been installed
## if not, install them
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

##############################################################
 #RUNNING STRUCTURE & Evanno Method
#Structure is meant to find the best K-Value for the population structure
#Helps to get most appropriate Q-Matrix
 structure <- function(strucFile, K.Max){
  
  strucFile <- toString(strucFile)
  print(strucFile)
  
  
  for(k in 1:K.Max){
      K <- toString(k)
      for(n in 1:10){
          x <- n
          x <- x + 10*(k-1)
          N <- toString(x)
          outfile <- paste("./StructureResults/results", "run",N, sep="_")
          outfile <- toString(outfile)
          system2("./structure", args = c("-i",strucFile,"-m","./mainparams", "-e", "./extraparams", "-K",K,"-o",outfile), wait=T)
      }
  }
  
  #RUNNING EVANNO METHOD TO DETERMINE APPROPRIATE K-VALUE
  sfiles <- list.files("./StructureResults", full.names = T)
  sfiles2 <- mixedsort(sfiles)
  slist <- readQ(sfiles2, filetype = "structure")
  tfiles <- tabulateQ(slist)
  sum <- summariseQ(tfiles)
  x <- evannoMethodStructure(sum,returnplot=T,returndata=F)
  
  jpeg("Evanno_Results.jpg")
  plot(x)
  dev.off()
  
  X11()
  plot(x)
  message("Press Return To Continue")
  invisible(readLines("stdin", n=1))
  cat("Enter chosen K-Value:")
  K <- readLines(file("stdin"), n=1)
  K.Val <- as.integer(K)
  K.Val <- K.Val*10
  df <- slist[K.Val]
  Qmatrix = data.frame(df)
  setwd("./Data")
  View(Qmatrix)
  write.table(Qmatrix, "./QMatrix_tmp.txt",quote=F,row.names=F)
}

##########################################################
#RUNNING GAPIT
#Both MLM and GLM methods are run and output is sent to GAPIT folder

gapit <- function(fileG, fileY, fileKI, K.Val){
  #GAPIT
  #Edited gapit source code, must have in current working directory
  setwd("../")
  source("./gapit_functions_mod.txt")
  
  #Genotype File (not numeric)  
  myG <- read.table(fileG, head= FALSE)
  
  #Phenotype File
  myY  <-read.table(fileY, head = TRUE)
  myY<-na_if(myY,"NH")
  write.table(myY, "./Data/Pheno.txt",quote=F,row.names=F)
  
  n <- ncol(myY)
  for (i in 2:n){
      myY[,i] <- as.numeric(myY[,i])
  }
  
  #K-Matrix
  myKI <- read.table(fileKI,head=FALSE)
  
  #Q-Matrix
  myCV <- read.table("./Data/QMatrix_tmp.txt",head=T)
  View(myCV)
  #myCV <- myCV[-nrow(myCV),]
  n <- ncol(myCV)
  n <- n+1
  myCV <- cbind(myY$ID, myCV)
  colnames(myCV)[1] <- 'TAXA'
  for (i in 2:n){
      val <- i-1
      colnames(myCV)[i] <- paste('Q',val,sep="")
  }
  #Need to drop last column, as GAPIT can fill in the blanks
  myC <- myCV[-ncol(myCV)]
  View(myC)
  
  #Running GAPIT  
  setwd("./GAPIT")
  View(myC)
  myGAPIT <- GAPIT(Y=myY,G=myG,KI=myKI, CV = myC, model=c("GLM","MLM"))
  
  View(myCV)
  colnames(myCV)[1] <- "<Trait>"
  write.table(myCV,"../Data/QMatrix.txt",quote=F,row.names=F)
}
########################################################
#Running tassel
#Runs both GLM & MLM and output is sent to TASSEL folder

tassel <- function(geno, phenoTassel, kinTassel){
#Running TASSEL
  setwd("../TASSEL")
  
   
  system("(echo \"<Covariate>\" && cat ../Data/QMatrix.txt) > ../Data/QMatrix_1.txt  && mv ../Data/QMatrix_1.txt  ../Data/QMatrix.txt",wait=T)

  
  #Running TASSEL GLM
  system2("../run_pipeline.pl", args = c("-fork1","-h", geno,"-filterAlign","-filterAlignMinCount","150", "-filterAlignMinFreq", "0.05", "-fork2","-r",phenoTassel, "-fork3","-q", "../Data/QMatrix.txt", "-combine4", "-input1","-input2","-input3","-intersect","-fork5","-k",kinTassel,"-combine6","-input5","-input4","-FixedEffectLMPlugin","-saveToFile","true", "-siteFile","tassel_stats_glm.txt","-alleleFile","tassel_allele_glm.txt","-endPlugin"), wait=T)
  
  #Running TASSEL MLM
  system2("../run_pipeline.pl", args = c("-fork1","-h", geno,"-filterAlign","-filterAlignMinCount","150", "-filterAlignMinFreq", "0.05", "-fork2","-r",phenoTassel, "-fork3","-q", "../Data/QMatrix.txt", "-combine4", "-input1","-input2","-input3","-intersect","-fork5","-k",kinTassel,"-combine6","-input4","-input3","-input5","-mlm","-mlmVarCompEst","P3D","-mlmOutputFile","tassel_results_mlm"), wait=T)
}
#######################################################

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
    df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    #select(chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    labs(y = "-log10(Q)") +
    
    # add genome-wide sig and sugg lines
    #geom_hline(yintercept = -log10(sig)) +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


#########################################################
#Merges results from both GAPIT & TASSEL
#Uses FDR method to determine most significant SNP's
#Merges results and plots new Manhattan Plots per stat method per phenotype
#Results can be found in Merged_Results folder


results <- function(columnNames){
   
  n <- length(columnNames)
  
  #setwd("../Merged_Results")
  setwd("./Merged_Results")

  #GAPIT & TASSEL RESULTS 
  for (i in 1:2){
    for(j in 2:n){
      pheno <- columnNames[j]
      gapit_name = ''
      tassel_name = ''
      tassel_df = ''
      if(i == 1){
        gapit_name <- paste("GAPIT.GLM.",pheno,".*\\.GWAS.Results.csv",sep='')
        tassel_name <- read.table("../TASSEL/tassel_stats_glm.txt",header=T)
        tassel_df <- data.frame(tassel_name)
      } else {
        gapit_name <- paste("GAPIT.MLM.",pheno,".*\\.GWAS.Results.csv",sep='')
        tassel_name <- list.files(path="../TASSEL/",pattern = ".*\\_stats.txt")
        tassel_file <- paste("../TASSEL/",tassel_name,sep='')
        tassel <- read.table(tassel_file, header=T, fill=T)
        tassel_df <- data.frame(tassel)
      }
      
      csv_file <- list.files(path="../GAPIT/",pattern = gapit_name)
      csv_file <- paste("../GAPIT/",csv_file,sep='')
      file <- read.csv(csv_file, head=T)
      gapit_df <- file[order(file[9]),]
      #Just include the SNP, Chromosome, Position, FDR P-Value
      gapit_df_shortened <- gapit_df[,c(1,2,3,9)]
      
      colnames(gapit_df_shortened)[4] <- "P"
      
      #Cut Off all Q-Values lower than 0.05
      gapit_cutOff <- gapit_df_shortened[which(gapit_df_shortened[4] <= 0.05),]

      
      #Reading in TASSEL results
      tassel_df <- tassel_df[,c(1,2,3,4,6)]
      colnames(tassel_df)[1] <- "Trait"
      colnames(tassel_df)[2] <- "SNP"
      colnames(tassel_df)[3] <- "CHR"
      colnames(tassel_df)[4] <- "BP"
      colnames(tassel_df)[5] <- "P"
      FDR_t = p.adjust(tassel_df$P,method="BH")
      tassel_df$FDR = FDR_t   
      tassel_df <- tassel_df[,-5]
      colnames(tassel_df)[5] <- "P"
      tassel_df_pheno <- tassel_df[which(tassel_df$Trait == pheno),]
      tassel_df_pheno <- tassel_df_pheno[order(tassel_df_pheno[5]),]
      tassel_df_cutOff <- tassel_df_pheno[which(tassel_df_pheno$P <= 0.05),]
      
      
      
      #Merging Results
      if(nrow(gapit_cutOff) > 0 & nrow(tassel_df_cutOff) > 0){
          merged = merge(x=tassel_df_cutOff, y=gapit_cutOff, by.x="SNP",by.y="SNP")
          #print(nrow(merged))
      
        if(nrow(merged) > 1){
            if(i == 1){
                fileName <- paste("Merged_GLM_",pheno,"_Stats.csv", sep='')
                plotName <- paste("ManhattanPlot_",pheno,"_GLM.jpg",sep='')
                name <-paste("ManhattanPlot_",pheno,"_GLM",sep='')
                ggplotName <- paste("ManhattanPlot_GG_",pheno,"_GLM.jpg",sep='')
                qqplotName <- paste("QQPlot_",pheno,"_GLM.jpg",sep='')
            } else{
                fileName <- paste("Merged_MLM_",pheno,"_Stats.csv", sep='')
                plotName <- paste("ManhattanPlot_",pheno,"_MLM.jpg",sep='')
                name <-paste("ManhattanPlot_",pheno,"_GLM",sep='')
                ggplotName <- paste("ManhattanPlot_GG_",pheno,"_GLM.jpg",sep='')
                qqplotName <- paste("QQPlot_",pheno,"_MLM.jpg",sep='')
            }
            #statresults
            write.csv(merged,fileName)
            
            #qq plot
            jpeg(qqplotName)
            qq(merged$P.y)
            dev.off()
       
       
            #Manhattan Plot
            mypalette <- c("#71D1CC", "#CB4577", "#0AA398", "#75002B", "#067E79") # chr color palette

            sig = 5e-8 # significant threshold line
            sugg = 1e-6 # suggestive threshold line
            ggMerged = merged[,c(1,3,4,5)]
            View(ggMerged)
            
            jpeg(ggplotName)
            colnames(ggMerged)[4] <- "P"
            print(gg.manhattan(ggMerged,threshold=1e-3,hlight=NA, ylim= c(0,5), title=name,col=mypalette))
            dev.off()
            
            mylabs <- unique(merged$CHR) # we need this to preserve the order, otherwise it will order alphabetically
            merged$CHR <- as.numeric(factor(merged$CHR, levels = mylabs))
            # then use custom labels
            #jpeg(plotName)
            #manhattan(merged,chr="CHR",chrlabs = merged$CHR,bp="Position",p="P.y",snp="SNP",logp=T,suggestiveline=F,genomewideline=F,main = "Manhattan Plot GLM", ylim = c(0, 5), ylab="-log(Q-Value)",cex = 0.6, cex.axis = 0.9, col = c("blue", "orange","red","purple","green"))
            #dev.off()
        } else {
            output=''
            if(i ==1)
                output <- paste("No Significant SNP's Found in ", pheno, " GLM ",  ,"Merged File!",sep='')
            else
                output <- paste("No Significant SNP's Found in ", pheno, " MLM ",  ,"Merged File!",sep='')
            print(output)
        }
      }
    }
  }
}

##########################################################
#Utilizing GLMNET for this and the lasso method,
#Binomial so phenotype is converted to discrete datatypes
#Uses Cross-Validation Method with 10 folds
#Plots are outputted to Lasso folder
#Coefficient results and such are printed to terminal


lasso <- function(columnNames,geno,myY){
    #Isolating SNPs data from Geno file
    setwd("../LASSO")
    geno <- read.table(geno,head=T)
    genoT = setNames(data.frame(t(geno[,-1])), geno[,1])
    genoTransposed <- cbind(names = rownames(genoT), genoT)
    rownames(genoTransposed) <- c()
    colnames(genoTransposed)[1] <- "SNP"
    
    
    n <- length(columnNames)
    
    #Need to go Through both GLM/MLM (outer loop)
    #Then inner loop is left for the phenotypes
    for (i in 1:2){
        for (j in 2:n){
            pheno <- columnNames[j]
            fileName =''
            if(i==1)
                fileName <- paste("../Merged_Results/Merged_GLM_",pheno,"_Stats.csv", sep='')
            else
                fileName <- paste("../Merged_Results/Merged_MLM_",pheno,"_Stats.csv", sep='')
            
            x <- -1
            if(file.exists(fileName)){
                currentPhenoResults <- read.csv(fileName)
                x <- nrow(currentPhenoResults)
            }
            
            
            #Check to Make sure there are SNP's in the read in file
            if(x > 1){
                genoJoin <- semi_join(genoTransposed, currentPhenoResults,by="SNP")
                genoJoinT = setNames(data.frame(t(genoJoin[,-1])), genoJoin[,1])
                genoJoinTransposed <- cbind(names = rownames(genoJoinT), genoJoinT)
                rownames(genoJoinTransposed) <- c()
                
                genoJoinTransposedM <- genoJoinTransposed[,-1]
                genoJoinMatrix <- as.matrix(genoJoinTransposedM)
                #genoJoinMatrix <- as.matrix(genoJoinTransposed)
                
                #Transforms phenotype into binary from numeric values, greater than 50% == 1
                currentPheno <- as.matrix(ifelse(myY[j]>50.0, 1, 0))

                
                #Running Lasso using Cross-Validation
                cv = cv.glmnet(genoJoinMatrix,currentPheno, family='binomial',nfolds=10,keep=T,trace.it=1,type.measure="class")
                print(coef(cv$glmnet.fit))
                print(cv)
                
                
                #fileCoef=''
                if(i==1){
                    fileCoeff <- paste(pheno,"_GLM_Coefficients_CV.txt",sep='')
                    fileJPG <- paste(pheno,"_GLM_CV_GLMNET_Coeffs.jpg",sep='')
                    fileJPG2 <- paste(pheno,"_GLM_CV_GLMNET.jpg",sep='')
                } else {
                    fileCoeff <- paste(pheno,"_MLM_Coefficients_CV.txt",sep='')
                    fileJPG <- paste(pheno,"_GLM_CV_GLMNET_Coeffs.jpg",sep='')
                    fileJPG2 <- paste(pheno,"_GLM_CV_GLMNET.jpg",sep='')
                }
                
                #Writing Files out to the Lasso Folder
                #cf <- as.data.frame(coef(cv$glmnet.fit))
                # cf <- tidy(coef(cv$glmnet.fit))
                #write.table(cf,fileCoeff)
                
                #Plotting Coefficients
                jpeg(fileJPG)
                plot(cv$glmnet.fit,label=T)
                dev.off()
                
                #Plotting Cross-Validation Plot
                jpeg(fileJPG2)
                plot(cv,label=T)
                dev.off()
            }
        }
    }
}

######################### MAIN ##################################

#format is only with one line per species, ploidy not specified

main <- function() {
  
  cat("\n\n---------------------------------------")
  cat("\n-------Running GWAS Pipeline-----------\n")
  cat("\n---------------------------------------\n\n")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  }  
  
  #### Install reqired packages in R Console
  packages<-c('readr','ggrepel','pophelper','LEA', 'fields', 'RColorBrewer', 'mapplots', 'dplyr', 'GAPIT3', 'gtools','devtools','ggplot2','gridExtra','gtable','label.switching','tidyr','qqman','glmnet')
  cat("\n\n*********Checking Dependency *********\n\n")
  check.packages(packages)

  
  cat("Enter Max K-Value to Test (will do an interval from 1:n):")
  K.Max <- readLines(file("stdin"), n=1)
  K.Max <- as.integer(K.Max)

  #Genotype file not numeric
  strucFile <- args[1]
  cat(strucFile)
  #structureFile <- read.table(strucFile,head=FALSE)
  
  #Running Structure
  cat("\n\nRunning Structure")
  #Qmatrix <- structure(strucFile, K.Max)
   
  
  #Running GAPIT
  fileG <- args[2]
  fileY <- args[3]
  fileKI <- args[4]
  cat("\n\nRunning GAPIT")
  #gapit(fileG,fileY,fileKI, K.Max)
  
  
  #Tassel Run
  geno <- paste(".",toString(fileG),sep="")
  phenoTassel <- paste(".",toString(args[5]),sep="")
  kinTassel <- paste(".",toString(args[6]),sep="")
  cat("\n\nRunning TASSEL")
  #tassel(geno, phenoTassel, kinTassel)
  
  #Number of Phenotypes
  getwd()
  #fileY <- paste(".",fileY,sep="")
  myY  <-read.table(fileY, head = TRUE)
  columnNames <- colnames(myY)
  
  #Results
  results(columnNames)
  cat("\n\nTASSEL and GAPIT Results have been merged and results can be found in Merged_Results Folder")
  
  #LASSO/GLMNET
  cat("\n\nRunning Lasso")
  numericGenoWithSNPNames <- paste(".",toString(args[7]),sep="")
  #lasso(columnNames,numericGenoWithSNPNames,myY)
  cat("\n\nLasso Results can be found in the Lasso Folder")
  
  #Done
  cat("\n\nPipeline Complete")
  cat("\n:)\n")
}

##################################################################


main()
