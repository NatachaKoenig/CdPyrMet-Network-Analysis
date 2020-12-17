#####################################################################
## PROJECT : ProteoGam                                             ##
## STUDIES : Tox network analysis of raw proteomic data            ##
## AUTHOR : Natacha Koenig                                         ##
## DATE : December 2020                                            ##
## SCRIPT : Network analysis of proteomics data from G. fossarum   ##
#####################################################################
#-------------------------------------------------------------------
#  INTRODUCTORY NOTE                            
#-------------------------------------------------------------------
# This script reproduce all the analytical steps used to produce 
# the results on the network analysis on the reproductive tissues of 
# Gammarus fossarum related to Cadmium - Pyriproxifen - Methoxyfenoside 
# contaminants.


# The sessionInfo() is in the appendice part

#-------------------------------------------------------------------
#  INDEX.SUMMARY OF THE PIPELINE                   
#-------------------------------------------------------------------

#      PACKAGES & SETTINGS 
#      DIRECTORIES  
#      VARIABLES DEFINTION
#      SCRIPT LOADING   

# 01.  DATA IMPORTATION
# 02.  FILTERING LOW ABUNDANCE PROTEINS OR NAs AND PRE-PROCESSING
# 03.  DATA NORMALIZATION 
# 04.  EXPLORATORY DATA ANALYSIS
# 05.  NETWORK CONSTRUCTION AND MODULE DETECTION
# 06.  RELATING MODULES TO TRAITS
# 07.  HUB GENES IN THE MODULES
# 08.  ANNOTATION
# 09.  VISUALIZATION OF THE PROTEIN NETWORK
# 10.  EXPORTING TO CYTOSCAPE  

#      APPENDICES

#-------------------------------------------------------------------
#  PACKAGES & SETTINGS                          
#-------------------------------------------------------------------
## Required packages
# Installs missing libraries !
list.of.packages <- c("plyr", "dplyr", "ggplot2", "grid", "gridExtra", "mixOmics", "minfi", "lumi", "stats", "WGCNA", "limma", "edgeR", "Heatplus", "made4", "RColorBrewer", "Biobase", "stringr") #list of packages required
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] #list of packages non installed
if(length(new.packages)) install.packages(new.packages, repos='https://cran.rstudio.com/') #install packages if the new.packages list is not empty

#tools
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

#PCA and clustering
library(mixOmics)
library(minfi)
library(dplyr)
library(lumi)
library(stats)

# Networks
library(WGCNA)
#If bug on blockwise function (namespace conflict)
cor <- WGCNA::cor
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Generics
library(limma)
library(edgeR)
library(Heatplus)
library(made4)
library(RColorBrewer)
library(Biobase)
library(stringr)


#-------------------------------------------------------------------
#  DIRECTORIES                          
#-------------------------------------------------------------------
# The script requires several sub-folders of the current directory :
# /data, /plot, /oitput_data and _cytoscape


## Working directory
wdir <- getwd()
wdir #current directory
dir()

## Input directories
datadir <- file.path(wdir, "data")

## Output directory
plotdir <- file.path(wdir, "plot")
output_datadir <- file.path(wdir, "output_data")
cytoscapedir <- file.path(output_datadir,"cytoscape")


#-------------------------------------------------------------------
#  VARIABLES DEFINITION                         
#-------------------------------------------------------------------
#To automate the analysis once adjusted
#we can define these variables upstream to launch the script at once

## PART 05. NETWORK CONSTRUCTION AND MODULE DETECTION 
# Definition of the newtwork construction parameters (l. 353)
softPower = 3 #soft-thresholding power for network construction
networkType="unsigned" #type of network
TOMType = "signed" #type of TOM
minModuleSize = 10 #the minimum Module Size 

## Choosing the method for the rest of the script and the filenames (l. 472)
# chosen_method = blockwise  #method of clustering #MUST BE DEFINED IN THE PART 5-NETWORK CONSTRUCTION AND MODULE DETECTION
method_for_file = "blockwise" #name of the method of clustering for the file extension


## PART 10. EXPORTING TO CYTOSCAPE
#All interactions are selected above this threshold (l. 739)
mod_threshold <- c(0.1, 0.05) #TOM threshold for cytoscape


#-------------------------------------------------------------------
#  SCRIPT LOADING                          
#-------------------------------------------------------------------
#the function script
source(file = "Functions_script.R")

#-------------------------------------------------------------------
#  01. DATA IMPORTATION                           
#-------------------------------------------------------------------
# Import the spectral count raw data of the tox dataset
rawdata <- read.delim(file.path(datadir, "tox-sc.txt"),check.names=FALSE, stringsAsFactors=FALSE)

# Taking a look to the raw data
names(rawdata) #samples names
row.names(rawdata)
rownames(rawdata) <- rawdata$GFOSS_contig_ID #replace the number of row by the contig names
rawdata<-rawdata[,-1] #delete the first column
head(rawdata)

# Read the phenotype data (i.e.: contaminant, exposure, class, group, batch effects)
pData <- read.delim(paste(datadir, "pData2.txt", sep="/"))
colnames(pData)
row.names(pData)<-pData$Sample_name #replace the number of row by the sample names
pData <- pData[,-1] #delete the first column
head(pData)

# Add of a column "Detailed_Class" to pData table to add the level of Control and Acetone
Detailed_Class <- pData$Class #Retrieve the "Class" column to a vector
Detailed_Class[6:10] <- rep("Acetone", 5) #replace "Control" by "Acetone"
pData[, "Detailed_Class"] <- as.factor(Detailed_Class) #Add the column "Detailed_Clss" to pData
head(pData)


#-------------------------------------------------------------------
#  02. FILTERING LOW ABUNDANCE PROTEINS OR MISSING VALUES      
#                           AND PRE-PROCESSING                      
#-------------------------------------------------------------------
## PRE-FILTERING STEP
# Put the data into a DGEList object (edgeR) y
y <- DGEList(counts=rawdata[1:40], genes=rownames(rawdata))
dim(y) # this object is still keeping proteins with too many NAs. In the next step we will filter
head(y$samples) #containing sample informations
head(y$genes) #containing annotation informations for each protein
save(y, file=file.path(output_datadir, "02_tox_y.RData"))

## LOOK AT THE DATA
## Euclidian clustering
pdf(file=file.path(plotdir, "02_CAH_before_filtering.pdf"), width=10, height=10)
plotSampleRelation(y$counts, subset = 871, method = "cluster", labels=pData$group, cex=1, main="Hierarchical clustering based on 871 proteins raw data")
dev.off()

## PCA
ty<-as.data.frame(t(y$counts)) #transpose count dataframe y

# Exploring how many dimensions explain variability
pdf(file=file.path(plotdir, "02_PCA_axes_before_filtering.pdf"))
tune.pca(ty, ncomp=10, center=T, scale=F) #here 3 dimensions
dev.off()

pca = pca(ty, ncomp=3) #we choose 3 dimensions for the PCA

# plotting the PCA results
pdf(file=file.path(plotdir, "02_PCA_before_filtering.pdf"))
plotIndiv(pca,
          group = pData$Group ,
          ind.names = F,
          legend=T, legend.position="right", legend.title=NULL,
          size.xlabel = 15,
          size.ylabel = 15,
          size.legend = 15,
          size.axis = 15,
          size.title = 15,
          title= paste('PCA rawdata samples'),
          comp=c(1,2),
          ellipse=F,
          pch=16
          )
dev.off()

## FILTERING
# Define a threshold corresponding to a count of 3 spectral counts for filtering
sc <- 3 #threshold = 3 spectral count
th <- cpm(sc, mean(colSums(rawdata)))[1] #converting SC to cpm
round(th) #the rounded threshold here is 2358
nbsamples <- 4 #mimimum number of samples threshold
keep <- rowSums(cpm(rawdata) > th) >= nbsamples #filtering step
rawdataF <- rawdata[keep,] #remaining proteins after filtering with "th" and "nbsamples"
dim(rawdataF) #312 remaining proteins
filtered_proteins <- row.names(rawdataF) #remaining proteins list


# Creation of two files for the filtered proteins list of the network and the associated spectral count (sc)
write.csv(filtered_proteins, file=file.path(output_datadir, "02_filtered_proteins_list.csv"))
write.csv(rawdataF, file=file.path(output_datadir, "02_filtered_proteins_sc.csv"))


#-------------------------------------------------------------------
#  03. DATA NORMALIZATION                      
#-------------------------------------------------------------------
# Normalization by the calcNormfactors function of the edgeR package
# based on a TMM (Trimmed mean of M-values) normalization procedure
#(https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)

# Put the filtered data into a DGEList object yF
yF<-DGEList(counts=rawdataF,
            genes=row.names(rawdataF),
            group= pData$Group)
dim(yF)
head(yF$samples) #containing sample informations
save(yF, file=file.path(output_datadir, "03_tox_yF.RData"))

# Normalization (required for WGCNA) yN
yN <- calcNormFactors(yF) #normalization of the filtered counts 
dim(yN)
head(yN$samples) #containing sample informations
save(yN, file=file.path(output_datadir, "03_tox_yN.RData"))

# Create the data with counts multiplied by yN norm.factors : yN2
yN2 <- yN
n <- length(colnames(yN$counts)) #number of columns of the yN counts table
n #number of samples

for (i in seq (from=1, to=n)) {
  yN2$counts[,i] <- yN$counts[,i]*yN$samples$norm.factors[i]
}

dim(yN2)
head(yN2$genes)
data <- yN2$counts #normalized counts data yN2 for the rest of the analysis
#rownames(data)<-yN2$genes$genes
rownames(data)
head(data)
save(yN2, data, pData, file = file.path(output_datadir, "03_tox_yN2.RData"))


#-------------------------------------------------------------------
#  04.  EXPLORATORY DATA ANALYSIS                   
#-------------------------------------------------------------------
##  LOOK AT THE DATA AFTER FILTERING AND NORMALIZATION
## Euclidian clustering
pdf(file=file.path(plotdir, "04_CAH_after_filt_norm.pdf"), width=10, height=10)
plotSampleRelation(yN2$counts, subset = 312, method = "cluster",labels=pData$Group, cex=1, main="Hierarchical clustering based on 312 proteins \n filtered and normalized data") #no more outliers
dev.off()

## PCA
tyN2<-as.data.frame(t(yN2$counts)) #transpose count dataframe yN2

# Exploring how many dimensions explain variability
pdf(file=paste(plotdir, "04_PCA_axes_after_filt_norm.pdf", sep="/"))
tune.pca(tyN2, ncomp=10, center=T, scale=F) #here 2 dimensions
dev.off()

pca =pca(tyN2, ncomp=2) #we choose 2 dimensions for the PCA

# plotting the PCA results
pdf(file=paste(plotdir, "04_PCA_after_filt_norm.pdf", sep="/"))
plotIndiv(pca,
          group = pData$Group ,
          ind.names = F,
          legend=T, legend.position="right", legend.title=NULL,
          size.xlabel = 15,
          size.ylabel = 15,
          size.legend = 15,
          size.axis = 15,
          size.title = 15,
          title= paste('PCA filtered and normalized samples'),
          comp=c(1,2),
          ellipse=F,
          pch=16)
dev.off()


#-------------------------------------------------------------------
#  05. NETWORK CONSTRUCTION AND MODULE DETECTION           
#-------------------------------------------------------------------
## Build the network based on proteins expressed over the previously established threshold

## CREATE AN ESET OBJECT to use for choosing the n. of genes to study with WGCNA.
esetpData<-as.data.frame(pData)
esetpData<-AnnotatedDataFrame(esetpData) #converting pData to an AnnotatedDataFrame
sampleNames(esetpData)=row.names(pData) #pData row names as sample names

eset<-ExpressionSet(data,esetpData) #containing expression level
eset

save(eset, file = paste(output_datadir, "05_tox_eset.RData", sep="/"))

## TRANSPOSE the table with expression values. 
# Here we use the whole protein dataset after the filtering step we performed above
data <- t(exprs(eset))
dim(data)

## MODEL THE DATA to approach a scale-free topology network for picking a soft-threshold.
# Choose a set of soft-thresholding powers
powers = c(c(1:10))
powers

# Call the network topology analysis function 
# From now, we have to work with the transposed table for expression values.
sft = pickSoftThreshold(data, powerVector = powers, verbose = 6)

# Plot the results (manual saving)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file=file.path(plotdir, "05_Soft_threshold_network.pdf"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,6],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6],
     labels=powers, cex=cex1,col="red")
#Plot line at the value of the mean connectivity which corresponds
#to the cutoff set previously
bestMean<-(sft$fitIndices[,6])
bestMean
#[1] 51.09056278 12.95741610  4.20421994  1.56731613  0.65572417  0.29476344  0.13937934
#[8]  0.07143150  0.03689017  0.01997293
abline (h=4.20421994 , col="red", lty =2)

dev.off()

## ONE-STEP NETWORK CONSTRUTION AND MODULE DETECTION
############# If these variables are not define above (variables section) uncomment ############ 
# # Definition of the newtwork construction parameters
# softPower = 3
# networkType="unsigned"
# TOMType = "signed"
# minModuleSize = 20

# Block-wise using all data
net = blockwiseModules(data, 
                       power = softPower, 
                       maxBlockSize=500,
                       networkType=networkType, 
                       TOMType = TOMType, 
                       minModuleSize = minModuleSize,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, 
                       deepSplit = 2,
                       saveTOMFileBase = "tox_net",
                       verbose = 3)

net # take a look at what looks like
save(net, file=file.path(output_datadir, "05_tox_net.RData"))

table(net$colors) #number of protein in each module (module 0 = unassigned proteins)
table(net$unmergedColors)
head(net$MEs) #preview of the ME related to module
net$goodSamples
table(net$goodGenes)
net$dendrograms
net$TOMFiles
table(net$blocks)
net$MEsOK

# Convert labels to colors for plotting
blockwiseColors = labels2colors(net$colors)
table(blockwiseColors)
# blockwiseColors
# black      blue     brown     green      grey   magenta      pink    purple       red turquoise    yellow 
# 23        42        39        24        38        16        20        12        23        45        30 

# Save the data
moduleLabels = net$colors
blockwiseColors = labels2colors(moduleLabels)
MEs = net$MEs
proteinTree = net$dendrograms[[1]]
save(MEs, moduleLabels, blockwiseColors, proteinTree, data, pData, eset, esetpData,
     file = file.path(output_datadir, "05_tox_network-blockwise.RData"))

# Plot the dendrogram and the module colors underneath (1 block)
# pdf(file=paste(plotdir, "PCA_after_filt_norm.pdf", sep="/"))
sizeGrWindow(10,5)
plotDendroAndColors(proteinTree, blockwiseColors,
                    "blockwise",
                    main = "Protein dendrogram and module definition",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()

## Module definition using TOM measures
# Refine the protein dendrogram, retrieving highly connected intramodular hub proteins
# Calculate  and store the adjacencies in an object, using the soft thresholding power 3.
ADJ1=abs(cor(data, use ="p"))^3
# or ADJ1= adjacency(data, power = softPower)
dissADJ1=1-ADJ1 #turn adjacency into a measure of dissimilarity

# Transform adjacency into Topologiacal Overlap Matrix (TOM) and calculate the corresponding dissimilarity
TOM=TOMsimilarity(ADJ1)
dissTOM=1-TOM

#Save data
write.csv(ADJ1, file.path(output_datadir, "05_tox_Adj.csv"))
write.csv(TOM, file.path(output_datadir, "05_tox_TOM.csv"))
write.csv(dissADJ1, file.path(output_datadir, "05_tox_dissADJ1.csv"))
write.csv(dissTOM, file.path(output_datadir, "05_tox_dissTOM.csv"))

# Module definition using TOM dissimilarity hierarchical clustering
# Build the dendrogram
hierTOM = hclust(as.dist(dissTOM), method="average")

colorStatTOM = as.character(cutreeStaticColor(proteinTree, cutHeight=.99, minSize=20))
colorDynTOM = labels2colors (cutreeDynamic(proteinTree,method="tree"))
colorDynHybrTOM = labels2colors(cutreeDynamic(proteinTree, distM=dissTOM , cutHeight=0.99,
                                             deepSplit=2, pamRespectsDendro=FALSE))
# Now we plot the results
pdf(file=paste(plotdir, "05_dendro_delim_modules.pdf", sep="/"), width=10)
# sizeGrWindow(10,5)
plotDendroAndColors(proteinTree,
                    colors=data.frame(blockwiseColors,colorStatTOM,
                                      colorDynTOM, colorDynHybrTOM),
                    dendroLabels = FALSE,
                    main = "Protein dendrogram & module definitions")
dev.off()

# Compare different module outcoume
table(blockwiseColors) ## blockwiseColors seems to be the best delilitation module method for our dataset
table(colorStatTOM)
table(colorDynTOM) 
table(colorDynHybrTOM)

# blockwiseColors
# black      blue     brown     green      grey   magenta      pink    purple       red turquoise    yellow 
# 23        42        39        24        38        16        20        12        23        45        30 

# colorStatTOM
# grey turquoise 
# 33       279 

# colorDynTOM
# blue     brown     green      grey       red turquoise    yellow 
# 49        40        26        85        20        58        34 

# colorDynHybrTOM
# blue     brown      grey turquoise    yellow 
# 77        72         1       129        33
  

############# If these variables are not define above (variables section) uncomment ############ 
# ## Choosing the method for the rest of the script and the filenames
chosen_method = blockwiseColors
# method_for_file = "blockwise"


## Calculating modules eigengenes in the chosen TOM model
model_ME=moduleEigengenes(data, chosen_method)$eigengenes
model_ME
#compare with blockwise calculated eigengenes
MEs

# Colors in the chosen meth
colList <- unique(chosen_method)
modules = colList


# Save the network showing good sensitivity and specificity, and relatively stringency in detecting genes with high MM
save(model_ME, blockwiseColors, modules, proteinTree, data, pData, eset, file = file.path(output_datadir, paste("05_tox_network_", method_for_file, ".RData", sep="")))



# Module summary
table(chosen_method)

# Loop to extract the proteins from each module of the chosen network 
lapply(1:length(modules), function(x){
  prot_list <- colnames(data)[chosen_method == modules[x]]
  print(length(prot_list))
  write.csv(prot_list, file.path(output_datadir, paste("05_proteins.", modules[x], "_tox-", method_for_file,".csv", sep="")))
})



#-------------------------------------------------------------------
#  06. RELATING MODULES TO TRAITS                     
#-------------------------------------------------------------------
# Load network data saved previously
lnames = load(file = file.path(output_datadir,
                               paste("05_tox_network_", method_for_file, ".RData", sep="")))
lnames

# Define numbers of proteins and samples
nGenes = ncol(data)
nGenes
nSamples = nrow(data)
nSamples

#Delete columns "Class", "Group" and "Detailed_Class"
head(pData)
pData2<-pData[,-c(1,2,14)] 
dim(pData2)

# Associate MEs with color labels
model_ME=moduleEigengenes(data, chosen_method)$eigengenes #calculating module eigengenes
MEs = orderMEs(model_ME) #put close eigenvextors next to each pther
moduleTraitCor = cor(MEs, pData2) #correlation value between MEs and samples
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) #pvalue related to the correlation of module's ME and samples

# Graphical representation of the module-trait correlation table
# We color code each association by the correlation value:
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot with [] gradient of contaminant
pdf(file=file.path(plotdir, "06_Modules_trait_network.pdf"), width=9, pointsize=14)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(pData2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = 0.8,
               cex.lab.y = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#### LIGHT VERSION  
# Associate MEs with color labels : without [] gradient of contaminant
moduleTraitCor_light=moduleTraitCor[,-c(4,5,7,8,10,11)] #correlation value between MEs and samples
moduleTraitPvalue_light=moduleTraitPvalue[,-c(4,5,7,8,10,11)] #pvalue related to the correlation of module's ME and samples
xLabels=names(pData2)[-c(4,5,7,8,10,11)]

# We color code each association by the correlation value:
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor_light, 2), "\n(",
                   signif(moduleTraitPvalue_light, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor_light)


# Display the correlation values within a heatmap plot
pdf(file=paste(plotdir, "06_Modules_trait_network_light.pdf", sep="/"), width=8, pointsize=16)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor_light,
               xLabels = xLabels,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab.x = 1,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



#-------------------------------------------------------------------
#  07. HUB GENES IN THE MODULES                   
#-------------------------------------------------------------------
## Compute intramodular connectivity for each protein
Alldegrees1=intramodularConnectivity(ADJ1, chosen_method) #kTotal, kWithin, kOut, kDiff
head(Alldegrees1)
write.csv(Alldegrees1, file.path(output_datadir, paste("07_tox_intramodularConnectivity_", method_for_file, ".csv", sep="")))

## Compute the protein Module Membership for each protein
geneModuleMembership=as.data.frame(cor(data, model_ME, use = "p")) #computing the protein Membership
head(geneModuleMembership)
# modNames = substring(names(geneModuleMembership), 3) #Retrieving the conames and deleting the "ME" at the beginning
# names(geneModuleMembership) = paste("MM", modNames, sep="") #ajoute "MM" devant les noms des colonnes qui correspondent aux couleurs des modules
write.csv(geneModuleMembership, file.path(output_datadir, paste("07_tox_geneModuleMembership_", method_for_file, ".csv", sep="")))

## Compute the pvalue for correlation related to the protein Module Membership for each protein
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
head(MMPvalue)
write.csv(MMPvalue, file.path(output_datadir, paste("07_tox_MMPvalue_", method_for_file, ".csv", sep="")))

## MDS plots to show proteins in modules with hub proteins in the finger tips
cmd1=cmdscale(as.dist(dissTOM),2) #classical multidimensional scaling (MDS)
sizeGrWindow(7,6)
par(mfrow=c(1,1))

pdf(file=file.path(plotdir, "07_MDS_plot_modules.pdf"))
plot(cmd1, col=as.character(chosen_method), main="MDS plot of network modules",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
dev.off()


#-------------------------------------------------------------------
#  08. ANNOTATION                
#-------------------------------------------------------------------
## DATA IMPORTATION
# annotation table
annotation_NCBI=read.csv(file=file.path(datadir, "annotation_NCBI_tox_NK.csv"), header=T, sep=';')
annotation_NCBI=tbl_df(annotation_NCBI) #create a dataframe tbl to apply the filter function

# kval importation (kTotal, kWithin, kOut, kDiff)
kval=read.csv(file=file.path(output_datadir, paste("07_tox_intramodularConnectivity_", method_for_file, ".csv", sep="")), header = T, sep=",")
head(kval)
colnames(kval)[1]="Contig" #replace the first column name by "Contig"

# MEval importation (the ModuleMembership values)
MEval=read.csv(file=file.path(output_datadir, paste("07_tox_geneModuleMembership_", method_for_file, ".csv", sep="")), header = T, sep=",")
head(MEval)
colnames(MEval)[1]="Contig" #replace the first column name by "Contig"


## NETWORK ANNOTATION
#sorting the module colors
modules_sorted <- sort(modules) #Sort module color by alphabetical order
modules_sorted

#Retrieve the protein list module and add the annotation associated
network_annot <- lapply(1:length(modules_sorted), function(x){
  
  prot_list <- colnames(data)[chosen_method == modules_sorted[x]]
  #print(length(prot_list))
  
  #filtering the big annotation table with the protein list of the module
  annot_filt_mod <- filter(annotation_NCBI, annotation_NCBI$Contig %in% prot_list)
  
  #adding the column of protein infos to the annot_filt_mod table : module color, kval and MEval
  full_annot_mod <- ADD_PROT_INFOS(annot_module=annot_filt_mod, color=modules_sorted[x], kval, MEval) #call the ADD_PROT_INFOS function
  
  #writing the full protein annotation table for each module
  write.csv(full_annot_mod, file=file.path(output_datadir, paste("08_annot_", modules_sorted[x], ".csv", sep="")))
  #write.csv(prot_list, file.path(output_datadir, paste("proteins.", modules_sorted[x], "_tox-", method_for_file,".csv", sep=""))) #c.f line 457 (part 5)
  return(full_annot_mod)
})


# CREATION OF AN ANNOTATION DATABASE TABLE WITH ALL MODULE ANOTATION
#do.call execute the function rbind to the list network_annot to create a table with all the protein module annotated
BDD_annot_tox <- do.call("rbind", network_annot) 
#export the database annotation network
write.csv(BDD_annot_tox, file=file.path(output_datadir, paste("08_BDD_annot_tox_", method_for_file, ".csv", sep=""))) 



#-------------------------------------------------------------------
#  09. VISUALIZATION OF THE PROTEIN NETWORK                
#-------------------------------------------------------------------
## Visualizing the protein network
# Load network data saved in the fifth part (construction of the network)
lnames = load(file = file.path(output_datadir, paste("05_tox_network_", method_for_file, ".RData", sep="")))
lnames #contains the names of loaded variables.

# Definition of the number of proteins and samples
dim(data) #number of protein (lignes) and samples (columns)
nGenes = ncol(data) #number of proteins
nSamples = nrow(data) #number of samples

# Load the topological overlap (TOM) for each block
dissTOM = 1-TOMsimilarityFromExpr(data, power = softPower) #the power define above (part 05)

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^9

# Set diagonal to NA  also improves the clarity of the plot
diag(plotTOM) = NA

# Call the plot function (heatmap of the network)
sizeGrWindow(9,9)
pdf(file=file.path(plotdir, paste("09_Heatmap_network_", method_for_file, ".pdf", sep="")))
TOMplot(plotTOM, proteinTree, chosen_method,
        main = "Network heatmap plot")
dev.off()

## Visualizing Eigengene relationships
# Plot the relationships among the eigengenes of the different modules (dendrogram + hm)
# Isolate contaminants and controls from the pData
Cd = as.data.frame(pData2$Cd)
names(Cd) = "Cd"
Met = as.data.frame(pData2$Met)
names(Met) = "Met"
Pyr = as.data.frame(pData2$Pyr)
names(Pyr) = "Pyr"

Control = as.data.frame(pData2$Control)
names(Control) = "Control"
Acetone = as.data.frame(pData2$Acetone)
names(Acetone) = "Acetone"

# Add the contaminants and the controls to existing module eigengenes
MET = orderMEs(cbind(model_ME, Cd, Met, Pyr, Control, Acetone))
sizeGrWindow(10,10)
par(cex = 0.9)
pdf(file=file.path(plotdir, paste("09_Eigengen_adjacency_", method_for_file, ".pdf", sep="")))
plotEigengeneNetworks(MET, "Eigengene adjacency", marDendro = c(1,4,2,2), 
                      marHeatmap = c(5,6,1,2), cex.lab = 0.8, 
                      xLabelsAngle= 90)
dev.off()



#-------------------------------------------------------------------
#  10. EXPORTING TO CYTOSCAPE                    
#-------------------------------------------------------------------
## Read the database module annotation
annotation_tox=read.csv(file=file.path(output_datadir, paste("08_BDD_annot_tox_", method_for_file, ".csv", sep="")), header=T, sep=',', row.names=1)

row.names(annotation_tox) <- annotation_tox$Contig #replacce the raw names by the name of protein
annotation_tox <- annotation_tox[,-1]
colnames(annotation_tox)[6] <- "Protein.name" #change the column "Name" for an other label because of a bug with Cytoscape importation
head(annotation_tox)


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(data, power = softPower) #the power defined above

############# If this variable is not define above (variables section) uncomment ############ 
# #All interactions are selected above this threshold
# mod_threshold <- c(0.1, 0.05)

# Select modules (colors in the chosen meth)
#modules <- unique(chosen_method)
#modules_sorted <- sort(modules)
modules_sorted

# Select module probes
probes = colnames(data) #names of proteins
inModule = is.finite(match(chosen_method, modules_sorted)) #check if the colors are in the chosen_method
modProbes = probes[inModule] #keep the probes which are inModule

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

## Export the network into edge and node list files Cytoscape can read
for (th in mod_threshold) {
  
  mod_cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = file.path(cytoscapedir, paste("10_CytoscapeInput-edges-",th, "-", paste(modules_sorted, collapse="-"), "_", 
                                                                        method_for_file, ".txt", sep="")),
                               nodeFile = NULL,
                               weighted = TRUE,
                               threshold = th,
                               nodeNames = modProbes,
                               nodeAttr = chosen_method[inModule])
  
  nodeFile <- merge(x = mod_cyt$nodeData, y=annotation_tox, by.x="nodeName", by.y="row.names")
  
  write.table(nodeFile, file.path(cytoscapedir, paste("10_CytoscapeInput-nodes-", th, "-", paste(modules_sorted, collapse="-"),  "_", 
                                                      method_for_file, ".txt", sep="")), sep="\t", quote=F, row.names=F)
  
}


#-------------------------------------------------------------------
#  APPENDICES                      
#-------------------------------------------------------------------
## sessionInfo()

# R version 3.5.0 (2018-04-23)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18362)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252 LC_NUMERIC=C                  
# [5] LC_TIME=French_France.1252    
# 
# attached base packages:
#   [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.4.0               made4_1.56.0                scatterplot3d_0.3-41        gplots_3.0.1.1             
# [5] RColorBrewer_1.1-2          ade4_1.7-13                 Heatplus_2.28.0             edgeR_3.24.3               
# [9] limma_3.38.3                WGCNA_1.68                  fastcluster_1.1.25          dynamicTreeCut_1.63-1      
# [13] lumi_2.34.0                 dplyr_0.8.3                 minfi_1.28.4                bumphunter_1.24.5          
# [17] locfit_1.5-9.1              iterators_1.0.12            foreach_1.4.7               Biostrings_2.50.2          
# [21] XVector_0.22.0              SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6        
# [25] matrixStats_0.55.0          Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
# [29] IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0         mixOmics_6.6.2             
# [33] lattice_0.20-35             MASS_7.3-49                 gridExtra_2.3               ggplot2_3.2.1              
# [37] plyr_1.8.4                 
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.4          Hmisc_4.2-0              igraph_1.2.4.1           lazyeval_0.2.2           splines_3.5.0           
# [6] robust_0.4-18.1          digest_0.6.21            htmltools_0.3.6          GO.db_3.7.0              gdata_2.18.0            
# [11] checkmate_1.9.4          magrittr_1.5             memoise_1.1.0            fit.models_0.5-14        cluster_2.0.7-1         
# [16] doParallel_1.0.15        readr_1.3.1              annotate_1.60.1          rARPACK_0.11-0           askpass_1.1             
# [21] siggenes_1.56.0          prettyunits_1.0.2        colorspace_1.4-1         rrcov_1.4-7              blob_1.2.0              
# [26] xfun_0.9                 crayon_1.3.4             RCurl_1.95-4.12          genefilter_1.64.0        impute_1.56.0           
# [31] GEOquery_2.50.5          zeallot_0.1.0            survival_2.41-3          glue_1.3.1               registry_0.5-1          
# [36] gtable_0.3.0             zlibbioc_1.28.0          Rhdf5lib_1.4.3           DEoptimR_1.0-8           HDF5Array_1.10.1        
# [41] scales_1.0.0             mvtnorm_1.0-11           DBI_1.0.0                rngtools_1.4             bibtex_0.4.2            
# [46] Rcpp_1.0.2               htmlTable_1.13.2         xtable_1.8-4             progress_1.2.2           foreign_0.8-70          
# [51] bit_1.1-14               mclust_5.4.5             preprocessCore_1.44.0    Formula_1.2-3            htmlwidgets_1.3         
# [56] httr_1.4.1               acepack_1.4.1            pkgconfig_2.0.3          reshape_0.8.8            XML_3.98-1.20           
# [61] nnet_7.3-12              labeling_0.3             tidyselect_0.2.5         rlang_0.4.0              reshape2_1.4.3          
# [66] AnnotationDbi_1.44.0     munsell_0.5.0            tools_3.5.0              RSQLite_2.1.2            knitr_1.25              
# [71] bit64_0.9-7              robustbase_0.93-5        beanplot_1.2             caTools_1.17.1.2         methylumi_2.28.0        
# [76] purrr_0.3.2              nlme_3.1-137             doRNG_1.7.1              nor1mix_1.3-0            xml2_1.2.2              
# [81] biomaRt_2.38.0           compiler_3.5.0           rstudioapi_0.10          affyio_1.52.0            tibble_2.1.3            
# [86] pcaPP_1.9-73             stringi_1.4.3            GenomicFeatures_1.34.8   RSpectra_0.15-0          Matrix_1.2-14           
# [91] multtest_2.38.0          vctrs_0.2.0              pillar_1.4.2             lifecycle_0.1.0          BiocManager_1.30.4      
# [96] data.table_1.12.2        bitops_1.0-6             corpcor_1.6.9            rtracklayer_1.42.2       latticeExtra_0.6-28     
# [101] R6_2.4.0                 affy_1.60.0              KernSmooth_2.23-15       nleqslv_3.3.2            codetools_0.2-15        
# [106] gtools_3.8.1             assertthat_0.2.1         rhdf5_2.26.2             openssl_1.4.1            pkgmaker_0.27           
# [111] withr_2.1.2              GenomicAlignments_1.18.1 Rsamtools_1.34.1         GenomeInfoDbData_1.2.0   mgcv_1.8-23             
# [116] hms_0.5.1                rpart_4.1-13             quadprog_1.5-7           tidyr_1.0.0              base64_2.0              
# [121] DelayedMatrixStats_1.4.0 illuminaio_0.24.0        base64enc_0.1-3          ellipse_0.4.1           