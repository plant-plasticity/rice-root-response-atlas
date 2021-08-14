##################################################################################
############ ATAC THS limma voom analysis, MDS plots, volcanos #################
##################################################################################
## Adapted from https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/DEG-analysis-limma-voom/Scripts/interactionDE-KK-SUB-SL.R 

######################## User defined options ########################
######################################################################
# Main directory:
setwd("~/mydirectory/") #Full path of the working directory. 
#It must contain: 
# ** A directory named "Counts" with a delimited file for raw counts. Rows are genes and columns samples.
# ** A directory named "meta" with a metadata file with information about the samples. Ideally the number of rows in the metadata is the same as in the raw counts.
# ** A directory named Scripts with this script and the 'functions.R' script.
## start a log file
sink('Log.txt')
## Metadata options
metaFile <- "metadata.csv" #Name of metadata file in csv format.
doFilter <- F #F
whichFilter <- c("") #If there are libraries that need to be filtered out (Avoid removing columns manually from the raw counts file)
## Counts file name (with extension) 
countsFile <- "ATAC_Counts_THSs.txt" #Name of counts file in tab delimited format
shortName<-"THS_date" #Name to add 
## Filter genes with low expression using CPM counts.
filterByCPM <- F #T
CPMcutoff <- 1 #2
## pValue (default = 0.05 ) and absolute logFC (default = 2) to color genes on volcano plots
pValCut=0.05  #0.05  
logCut=2 #2
########################
########################

###
library(edgeR)
library(reshape)
library(gplots)
library(RColorBrewer)
library(calibrate)
library(Glimma) #source("https://bioconductor.org/biocLite.R")biocLite("Glimma")

## Output
outDir = "Outdir/"
dir.create(outDir, showWarnings=T)

geneListsDir = "Outdir/GeneLists"
dir.create(geneListsDir, showWarnings=T)
#
imgDir = "Outdir/images/"
dir.create(imgDir, showWarnings=T)
## --

if (is.na(shortName)){
  shortName <- basename(getwd())
}

# Load functions
source("Scripts/functions.R")
######## --- --- --- 


## Start of analysis
####################################################################################
####################################################################################
cat("Reading metadata file \n")

meta <- metaDataProcessing(metaFile,doFilter)
head(meta)

#
cat("Reading counts file:",countsFile,"\n")

GeneCounts <- read.delim(paste0("Counts/",countsFile),row.names = 1)
dim(GeneCounts)

## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
tmp <- filterCounts(GeneCounts,meta)
GeneCounts <- tmp[["counts"]]
meta <- tmp[["meta"]]
rm(tmp)
## --


###### Design matrix
## Convert experimental metadata to factors for the design
experimentFactors <- lapply(apply(meta,2,split,""),unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors,as.factor))

cat ("Create the design with these factors:\n")
print(head(experimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
Groups <- as.factor(paste0(experimentFactors$Sample,experimentFactors$Treatment,experimentFactors$Genotype,experimentFactors$Tissue))
design <- model.matrix(~0+Groups) 

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","experimentFactors","\\$","\\:","\\-",
                   colnames(experimentFactors)),sep="",collapse = "|")

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)


####################################################################################
cat("Removing THSs with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
cat("Removing",length(rmIDX),"genes \n")
GeneCounts <- GeneCounts[-rmIDX,]
cat("Remaining number of THS:",nrow(GeneCounts),"\n")


### Use cpms to uncover lowly expressed genes
dge <- DGEList(counts=GeneCounts,remove.zeros = T)

# Filter genes with low CPMs accross replicates 
cat("Replicates of samples range between:", range(table(Groups)),"\n")

#
if (filterByCPM){
  
  sampleMin <- min(table(Groups))
  cat("Filtering reads with low CPMs ( <",CPMcutoff,") in at least",sampleMin,"replicates \n")
  #
  cpm <- cpm(dge)
  keep.exprs <- rowSums(cpm>CPMcutoff)>=sampleMin
  table(keep.exprs)
  
  
  cat("Removing",table(keep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(keep.exprs)[2],"\n")
  
  #
  
  y <- dge[keep.exprs, , keep.lib.size = FALSE]
} else {
  cat("Not doing CPM filtering")
}

normalizedExpression <- cpm(y)
#
tmpSave <- paste(outDir,"normalizedExpression","_",shortName,".csv",sep="")
cat("Saving normalized data to:", tmpSave, "\n")
save(normalizedExpression,file = "cpm_normalizedExpression.RData")
write.csv(x=normalizedExpression,tmpSave,quote = F,row.names = T)


## Easier visualization for MDS plots
cat("Using glimma for MDS plot visualization - Normalized data \n")
glMDSPlot(y, labels=rownames(y$samples),
          groups=meta,folder = "ATAC_THSs/", launch=T)

#### Start PDF
tmpSave <- paste(imgDir,"DEG_Analysis_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")

### Use voom on the dge object.
v <- voom(y, design, plot = TRUE,normalize.method ="quantile")
# OR:
#v <- voomWithQualityWeights(y,design, normalization="quantile",plot = T)
###
cat("Analyzing",nrow(v),"with",ncol(v),"libraries \n")

## Obtain back quantile normalized reads
r=v
indsamp=length(colnames(r$E))
r$E[,1:indsamp]<-2^r$E[,1:indsamp] ##revert log
## calculate million reads
m<-r$targets$lib.size/1000000
## transform reads
r$E=t(t(r$E)*m)

## get Granges for THSs
Features <- readRDS("THSs.RDS")
## names
names(Features)=as.character(1:length(Features))
## Features object over remaining THSs
Features_present=Features[names(Features)%in%rownames(r$E)]

## RPKM calculation
library(systemPipeR)
rpkmDFeByg <- apply(r$E, 2, function(x) returnRPKM(counts=x, ranges=Features_present))
rpkmDFeByg=rpkmDFeByg[,order(colnames(rpkmDFeByg))]

## TPM calculation
sums=colSums(rpkmDFeByg)
TPM=rpkmDFeByg
for (i in 1:length(colnames(rpkmDFeByg))){
  TPM[,i]=rpkmDFeByg[,i]/sums[i]*10^6
}
## use mean function
meanRPKM <- meanNormalizedExpression(rpkmDFeByg,levels(Groups)) 
meanTPM <- meanNormalizedExpression(TPM,levels(Groups)) 
write.table(rpkmDFeByg, paste0(outDir,"RPKM_OS_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(TPM, paste0(outDir,"TPM_OS_after_voom_p_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(meanRPKM, paste0(outDir,"RPKM_OS_Mean_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(meanTPM, paste0(outDir,"TPM_OS_Mean_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
######## Visualization and quality control
#testPalette(Colors13)
#testPalette(ColoresPair)
#testPalette(customColors) #18 colors

##################
## Correlation between replicates of samples belonging to same group
corrSamples <- cor(v$E)
## --
tmpSave <- paste(imgDir,"CorrelationBetweenReplicates_p_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")#width = 8,height = 6)
#colors <- colorRampPalette(c("darkgoldenrod4","darkgoldenrod1","white","white","steelblue1","steelblue4"))

for (each in (levels(Groups))){
  hmData <- corrSamples[grep(each, rownames(corrSamples)),grep(each, colnames(corrSamples))]
  #hmData <- corrSamples[,grep(each, colnames(corrSamples))]
  hm <- T
  if(!hm){
    cat("Heatmaps with NMF \n")
    NMF::aheatmap(hmData,col=colors(125),
                  txt = ifelse(hmData<0.85,"<",NA),#Rowv = F,Colv = F,
                  main=paste0("Correlation between samples of group ",each))
    
  } else {
    cat("Heatmaps with heatmap.2 \n")
    heatmap.2(hmData,col=cm.colors(125), keysize = 0.75,
              cellnote = ifelse(hmData<0.8,"*",NA), notecol = "black",
              
              #margins = c(16,16),
              
              dendrogram = "none", trace = "none",density.info='none',
              cexCol  = 0.8 ,cexRow = 0.8,
              lmat=rbind(c(4, 3, 9),
                         c(2, 1, 6),
                         c(8, 5, 7)),
              lhei=c(0.3, 0.6,0.8),
              lwid=c(0.25, 0.4,0.2),
              main=paste0("Correlation\n",each))
    legend("bottomleft",legend = "* means correlation < 0.85",bty = "n")
  }
  
}
dev.off()

write.csv(corrSamples[order(rownames(corrSamples)),order(colnames(corrSamples))],"corrSamples.csv")


## Assign colors to each of the experimental factors. 
ColorTable <- assignColorsToFactors(experimentFactors)

## Boxplot of normalized counts ordered by Groups
boxplot(v$E[,order(Groups)], range=0,col=customColors[Groups[order(Groups)]], 
        ylab="log2[counts]", xlab="sample", main="Quantile normalized Counts",
        cex.axis=0.5,las=2)


## Contrast for comparisons
cont.matrix= makeContrasts(
  "ATSUB"=ATSB-ATND, #should be adjusted based on naming for groups 
  "ATSUBREC"=ATSUBREC-ATSB,
  "ATSUBREC_CON"=ATSUBREC-ATND,
  "ATWD"=ATDO-ATCON,
  "ATWDREC"=ATDROREC-ATDO,
  "ATWDREC_CON"=ATDROREC-ATCON,
  "ATWL"=ATWT-ATCON,
  "ATWLREC"=ATWATREC-ATWT,
  "ATWLREC_CON"=ATWATREC-ATCON,
  "ATFIELD"=ATFIELD-ATCON,
  "ATFIELD_CON"=ATFIELD-ATND,
  "ATCON_ND"=ATCON-ATND,
  "ATFIELD_WL"=ATFIELD-ATWT,
  "ATFIELD_WLREC"=ATFIELD-ATWATREC,
  "ATFIELD_WD"=ATFIELD-ATDO,
  "ATFIELD_WDREC"=ATFIELD-ATDROREC,
  "ATFIELD_SUBREC"=ATFIELD-ATSUBREC,
  "ATFIELD_SUB"=ATFIELD-ATSB,
  levels=design)

#### Fit and do differential accessibility calculation
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## -- Summary and Venn diagrams , only good for up to 5 comparisons.
results <- decideTests(fit2)
summary(results)
if (ncol(results) <= 5){
  cat ("Doing Venn Diagrams \n")
  vennDiagram(results,include = c("up","down"), main="DE")
} else {
  cat ("More than 5 comparisons, skipping Venn Diagrams  \n")
}
DESummary <- t(summary(decideTests(fit2)))[,-2]
colnames(DESummary) = c("Downregulated","Upregulated")

# Save as csv
tmpSave <- paste(outDir,"DESummaryInteractions_",shortName,".csv",sep="")
write.csv(x=DESummary,tmpSave,quote = F,row.names = T)

# Write to PDF
plotData <- t(DESummary)
yMax <- max(colSums(plotData))
rownames(plotData) <- c("Down","Up")

barplot(plotData,legend.text = rownames(plotData),col=c("orange","steelblue4"),
        xlab = "Contrast", ylab = "Number of genes",
        beside = T,
        ylim = c(0,yMax*1.2),
        las=2,
        cex.names = 0.6, border = T, bty="n",
        main="DE genes per contrast")

DEList <- list()
for (contrast in colnames(cont.matrix)){
  print(contrast)
  ## Sorting by none ensures all contrasts will be in the same order
  tmp <- topTable(fit2, coef=contrast,number = Inf,sort.by = "none")
  #
  pValpassed <- table(tmp$adj.P.Val < 0.05)[2]
  cat ("Number of genes with adj pVal < 0.05 on ",contrast,":",pValpassed,"\n")
  
  
  ## Write genes that are up or downregulated (logFC > 1; logFC < (-1))
  upGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC > 1,]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_up",".csv",sep="")
  write.csv(x=upGenes,tmpSave,quote = F,row.names = T)
  #
  downGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC < (-1),]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_down",".csv",sep="")
  write.csv(x=downGenes,tmpSave,quote = F,row.names = T)
  #####
  
  #-- Add gene symbols if available
  tmp[,"Symbol"] <- rownames(tmp)
  if (annotationAvail){
    cat("Adding annotation \n")
    Genes <- rownames(tmp)
    idx <- intersect(names(AGI2Symbol),Genes)
    tmp[idx,"Symbol"] <- AGI2Symbol[idx]
    Genes
  }
  #--
  
  ## Add contrast name to the column names, in case of multiple contrasts.
  colnames(tmp) <- paste(colnames(tmp),contrast,sep = ".")
  
  # Write each contrast to file
  tmpSave <- paste(outDir,contrast,"_",shortName,".csv",sep="")
  write.csv(x=tmp,tmpSave,quote = F,row.names = T)
  
  # Save result to list
  DEList[[contrast]] <- tmp 
}

tmpSave <- paste(outDir,"DEList_",shortName,".RData",sep="")
save(DEList,file = tmpSave)
### ------------


tmpSave <- paste(imgDir,"VolcanoPlots_0.05_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")
makeVolcanoPlots(DEList,pValCut=0.05,logCut=2,plotGenes=F) #plotGenes=T to print genes in the plot
dev.off()

### Condense into a single list
## Convert to a single table
DE_All <- condenseListTables(DEList) ## Use a custom function

DE_All <- DE_All[,-grep("t.|B.|P.Value|AveExpr",colnames(DE_All))] #Remove unwanted columns

saveRDS(DE_All,"All_contrasts.RDS")

THS_counts=read.csv("/mydirectory/THS_GH_counts_annotated.csv")
head(THS_counts)
THS_counts[as.numeric(names(eByg2)),30:44]
##
DE_All <- cbind(DE_All,meanTPM[rownames(DE_All),],THS_counts[as.numeric(names(eByg2)),30:44])

## Write output
tmpSave <- paste(outDir,"DEG_AllContrasts_",shortName,".csv",sep="")
write.csv(x = DE_All,file = tmpSave,quote = F,row.names = T)

## Close main img
dev.off()
dev.off()
sink()

######################################################################################
