##################################################################################
############ polyA / TRAP limma voom analysis, MDS plots, volcanos #################
##################################################################################
## Adapted from https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/DEG-analysis-limma-voom/Scripts/interactionDE-KK-SUB-SL.R 


######################## User defined options ########################
######################################################################
# Main directory:
setwd("/Users/mauricio/Documents/bioinformatics_meetings/limmavoom_v2_may10") #Full path of the working directory. 
#It must contain: 
# ** A directory named "Counts" with a delimited file for raw counts. Rows are genes and columns samples.
# ** A directory named "meta" with a metadata file with information about the samples. Ideally the number of rows in the metadata is the same as in the raw counts.
# ** A directory named Scripts with this script and the 'functions.R' script.

sink('Log_file.txt')
## Metadata options
metaFile <-"metadata.csv"  #Name of metadata file
doFilter <- F #F
whichFilter <- c("TOTCONHAIRYROOTT6A") #If there are libraries that need to be filtered out (Avoid removing columns manually from the raw counts file)


countsFile<-"ALL_Counts.txt" ## Counts file name (with extension) 

shortName <- "GH_Field_Plate_polyA_TRAP_date" #Short name to append at the end of the filenames. If missing it will append the name of the folder where the scripts where run.

## Filter genes with low expression using CPM counts.
filterByCPM <- T #T
CPMcutoff <- 3 #2


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
outDir = "GH_Field_Plate_polyA_TRAP_date/"
dir.create(outDir, showWarnings=T)

geneListsDir = paste0 (outDir,"GeneLists/")
dir.create(geneListsDir, showWarnings=T)
#
imgDir = paste0 (outDir,"images/")

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
Groups <- as.factor(paste0(experimentFactors$Promoter,experimentFactors$Sample,experimentFactors$Treatment,experimentFactors$Genotype,experimentFactors$Tissue))
# Groups
design <- model.matrix(~0+Groups) 
# Example of an interaction
#design <- model.matrix(~0+experimentFactors$Sample*experimentFactors$Treatment) #Sample*Treatment interaction

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","experimentFactors","\\$","\\:","\\-",
                   colnames(experimentFactors)),sep="",collapse = "|")

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)


####################################################################################
cat("Removing genes with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
cat("Removing",length(rmIDX),"genes \n")
GeneCounts <- GeneCounts[-rmIDX,]
cat("Remaining number of genes:",nrow(GeneCounts),"\n")


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

## CPM calculation
normalizedExpression <- cpm(y)
#
###### Adjustments to get with the same data the MDS plot only from 35S samples
## Easier visualization for MDS plots

cat("Using glimma for MDS plot visualization - Normalized data \n")
glMDSPlot(y, labels=rownames(y$samples),
          groups=meta,folder=paste0("glimma_",shortName), launch=T)


#### Start PDF
tmpSave <- paste(imgDir,"DEG_Analysis_",shortName,".pdf",sep="")

pdf(tmpSave,paper = "USr")


## Use voom on the dge object apply quantile normalization
v <- voom(y, design, plot = TRUE,normalize.method ="quantile")
# OR:
#v <- voomWithQualityWeights(y,design, normalization="quantile",plot = T)
###
cat("Analyzing",nrow(v),"with",ncol(v),"libraries \n")

## Obtain back quantile normalized reads
r=v
indsamp=length(colnames(r$E)) ##number of columns
r$E[,1:indsamp]<-2^r$E[,1:indsamp] ##revert log
## calculate million reads
m<-r$targets$lib.size/1000000
## transform reads
r$E=t(t(r$E)*m)

library(GenomicFeatures); library(systemPipeR)
txdb=loadDb("/txdb.sqlite") #transcript database object
eByg <- exonsBy(txdb, by="gene") ## get Granges for exonic regions
rownames(r$E)=substr(rownames(r$E),1,12) ## apply rownames
eByg_selection=eByg[names(eByg)%in%rownames(r$E)] ## extract only genes calculated
## calculate RPKM
rpkmDFeByg <- apply(r$E, 2, function(x) returnRPKM(counts=x, ranges=eByg_selection))
rpkmDFeByg=rpkmDFeByg[,order(colnames(rpkmDFeByg))]
## calculate TPM
sums=colSums(rpkmDFeByg)
TPM=rpkmDFeByg
for (i in 1:length(colnames(rpkmDFeByg))){
  TPM[,i]=rpkmDFeByg[,i]/sums[i]*10^6
}

meanRPKM <- meanNormalizedExpression(rpkmDFeByg,levels(Groups)) 
meanTPM <- meanNormalizedExpression(TPM,levels(Groups)) 
anno=readRDS("annotation_040620.RDS")
head(anno)
## function to add annotation
annotate = function (File_of_interest,anno){
  File_of_interest=as.data.frame(File_of_interest)
  File_of_interest[,(ncol(File_of_interest)+1):(ncol(File_of_interest)+ncol(anno)-1)]=anno[match(rownames(File_of_interest),substr(anno$GENEID,1,12)),2:ncol(anno)]
  return(File_of_interest)
}
rpkmDFeByg=annotate(rpkmDFeByg,anno)
TPM=annotate(TPM,anno)
meanRPKM=annotate(meanRPKM,anno)
meanTPM=annotate(meanTPM,anno)

## Export
write.table(rpkmDFeByg, paste0(outDir,"RPKM_OS_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(TPM, paste0(outDir,"TPM_OS_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(meanRPKM, paste0(outDir,"RPKM_OS_Mean_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")
write.table(meanTPM, paste0(outDir,"TPM_OS_Mean_after_voom_",shortName,".xls"), col.names=NA, quote=FALSE, sep="\t")

######## Visualization and quality control
#testPalette(Colors13)
#testPalette(ColoresPair)
#testPalette(customColors) #18 colors

##################
## Correlation between replicates of samples belonging to same group
corrSamples <- cor(v$E)

my_palette <- colorRampPalette(c("darkgoldenrod4","darkgoldenrod1","white","white","steelblue1","steelblue4"))
c2=corrSamples[,order(colnames(corrSamples))]
c2=c2[order(rownames(c2)),]
## --
tmpSave <- paste(imgDir,"CorrelationBetweenReplicates_7_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")#width = 8,height = 6)
#colors <- colorRampPalette(c("darkgoldenrod4","darkgoldenrod1","white","white","steelblue1","steelblue4"))

for (each in (levels(Groups))){
  hmData <- corrSamples[grep(each, rownames(corrSamples)),grep(each, colnames(corrSamples))]
  #hmData <- corrSamples[,grep(each, colnames(corrSamples))]
  hm <- T
  if(!hm){
    cat("Heatmaps with NMF \n")
    NMF::aheatmap(hmData,col=my_palette,
                  txt = ifelse(hmData<0.8,"<",NA),#Rowv = F,Colv = F,
                  main=paste0("Correlation between samples of group ",each))
    
  } else {
    cat("Heatmaps with heatmap.2 \n")
    heatmap.2(hmData,col=my_palette, keysize = 0.75,
              cellnote = ifelse(hmData<0.8,"*",NA), notecol = "black",
              
              #margins = c(16,16),
              breaks=seq(0.7,1,by=0.3/length(my_palette)),
              dendrogram = "none", trace = "none",density.info='none',
              cexCol  = 0.8 ,cexRow = 0.8,
              lmat=rbind(c(4, 3, 9),
                         c(2, 1, 6),
                         c(8, 5, 7)),
              lhei=c(0.3, 0.6,0.8),
              lwid=c(0.25, 0.4,0.2),
              main=paste0("Correlation\n",each))
    legend("bottomleft",legend = "* means correlation < 0.8",bty = "n")
  }
  
}
dev.off()

################## --
tmpSave <- paste(imgDir,"CorrelationBetweenReplicates_all_",shortName,".pdf",sep="")
pdf(tmpSave,width = 18,height = 18)
if(!hm){
  cat("Heatmaps with NMF \n")
  NMF::aheatmap(c2,col=my_palette,
                txt = ifelse(hmData<0.8,"<",NA),#Rowv = F,Colv = F,
                main=paste0("Correlation between samples of group ",each))
  
} else {
  cat("Heatmaps with heatmap.2 \n")
  heatmap.2(c2,col=my_palette, keysize = 0.2,
            cellnote = ifelse(c2<0.8,"*",NA), notecol = "black",
            
            #margins = c(16,16),
            breaks=seq(0.7,1,by=0.3/length(my_palette)),
            dendrogram = "none", trace = "none",density.info='none',
            cexCol  = 0.8 ,cexRow = 0.8,
            lmat=rbind(c(4, 3, 9),
                       c(2, 1, 6),
                       c(8, 5, 7)),
            lhei=c(0.3, 0.6,0.8),
            lwid=c(0.25, 0.4,0.2),
            main=paste0("Correlation\n",each))
  legend("bottomleft",legend = "* means correlation < 0.8",bty = "n")
}
dev.off()


tmpSave <- paste(imgDir,"Extra_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")
## Assign colors to each of the experimental factors. 
ColorTable <- assignColorsToFactors(experimentFactors)

## Boxplot of normalized counts ordered by Groups
boxplot(v$E[,order(Groups)], range=0,col=customColors[Groups[order(Groups)]], 
        ylab="log2[counts]", xlab="sample", main="Quantile normalized Counts",
        cex.axis=0.5,las=2)


### Contrast tables calculated as groups

# group 1
cont.matrix= makeContrasts(
  "X35SWD-CON-TRAP"=X35STRADROOSROOT-X35STRACONOSROOT,
  "X35SWD1R-CON-TRAP"=X35STRADRORECOSROOT-X35STRACONOSROOT,
  "X35SWD1R-WD-TRAP"=X35STRADRORECOSROOT-X35STRADROOSROOT,
  "X35SWL-CON-TRAP"=X35STRAWATOSROOT-X35STRACONOSROOT,
  "X35SWL1R-CON-TRAP"=X35STRAWATRECOSROOT-X35STRACONOSROOT,
  "X35SWL1R-WL-TRAP"=X35STRAWATRECOSROOT-X35STRAWATOSROOT,
  "X35SSUB-CON-TRAP"=X35STRASUBOSROOT-X35STRACON2OSROOT,
  "X35SSUB1R-CON-TRAP"=X35STRASUBRECOSROOT-X35STRACON2OSROOT,
  "X35SSUB1R-SUB-TRAP"=X35STRASUBRECOSROOT-X35STRASUBOSROOT,
  "X35SFIELD-CON-TRAP"=X35STRAFIELDOSROOT-X35STRACONOSROOT,
  "X35SPLATE-CON-TRAP"=X35STRAPLATEOSROOT-X35STRACONOSROOT,
  "LSI1WD-CON-TRAP"=LSI1TRADROOSROOT-LSI1TRACONOSROOT,
  "LSI1WD1R-CON-TRAP"=LSI1TRADRORECOSROOT-LSI1TRACONOSROOT,
  "LSI1WD1R-WD-TRAP"=LSI1TRADRORECOSROOT-LSI1TRADROOSROOT,
  "LSI1WL-CON-TRAP"=LSI1TRAWATOSROOT-LSI1TRACONOSROOT,
  "LSI1WL1R-CON-TRAP"=LSI1TRAWATRECOSROOT-LSI1TRACONOSROOT,
  "LSI1WL1R-WL-TRAP"=LSI1TRAWATRECOSROOT-LSI1TRAWATOSROOT,
  "LSI1FIELD-CON-TRAP"=LSI1TRAFIELDOSROOT-LSI1TRACONOSROOT,
  "RSS1WD-CON-TRAP"=RSS1TRADROOSROOT-RSS1TRACONOSROOT,
  "RSS1WD1R-CON-TRAP"=RSS1TRADRORECOSROOT-RSS1TRACONOSROOT,
  "RSS1WD1R-WD-TRAP"=RSS1TRADRORECOSROOT-RSS1TRADROOSROOT,
  "RSS1WL-CON-TRAP"=RSS1TRAWATOSROOT-RSS1TRACONOSROOT,
  "RSS1WL1R-CON-TRAP"=RSS1TRAWATRECOSROOT-RSS1TRACONOSROOT,
  "RSS1WL1R-WL-TRAP"=RSS1TRAWATRECOSROOT-RSS1TRAWATOSROOT,
  "RSS1SUB-CON-TRAP"=RSS1TRASUBOSROOT-RSS1TRACON2OSROOT,
  "RSS1SUB1R-CON-TRAP"=RSS1TRASUBRECOSROOT-RSS1TRACON2OSROOT,
  "RSS1SUB1R-SUB-TRAP"=RSS1TRASUBRECOSROOT-RSS1TRASUBOSROOT,
  "RSS1FIELD-CON-TRAP"=RSS1TRAFIELDOSROOT-RSS1TRACONOSROOT,
  "RSS1PLATE-CON-TRAP"=RSS1TRAPLATEOSROOT-RSS1TRACONOSROOT,
  "QHBWD-CON-TRAP"=QHBTRADROOSROOT-QHBTRACONOSROOT,
  "QHBFIELD-CON-TRAP"=QHBTRAFIELDOSROOT-QHBTRACONOSROOT,
  "QHBPLATE-CON-TRAP"=QHBTRAPLATEOSROOT-QHBTRACONOSROOT,
  "CMZWD-CON-TRAP"=CMZTRADROOSROOT-CMZTRACONOSROOT,
  "CMZWD1R-CON-TRAP"=CMZTRADRORECOSROOT-CMZTRACONOSROOT,
  "CMZWD1R-WD-TRAP"=CMZTRADRORECOSROOT-CMZTRADROOSROOT,
  "CMZWL-CON-TRAP"=CMZTRAWATOSROOT-CMZTRACONOSROOT,
  "CMZWL1R-CON-TRAP"=CMZTRAWATRECOSROOT-CMZTRACONOSROOT,
  "CMZWL1R-WL-TRAP"=CMZTRAWATRECOSROOT-CMZTRAWATOSROOT,
  "CMZPLATE-CON-TRAP"=CMZTRAPLATEOSROOT-CMZTRACONOSROOT,
  "CASPWD-CON-TRAP"=CASPTRADROOSROOT-CASPTRACONOSROOT,
  "CASPWD1R-CON-TRAP"=CASPTRADRORECOSROOT-CASPTRACONOSROOT,
  "CASPWD1R-WD-TRAP"=CASPTRADRORECOSROOT-CASPTRADROOSROOT,
  "CASPWL-CON-TRAP"=CASPTRAWATOSROOT-CASPTRACONOSROOT,
  "CASPWL1R-CON-TRAP"=CASPTRAWATRECOSROOT-CASPTRACONOSROOT,
  "CASPWL1R-WL-TRAP"=CASPTRAWATRECOSROOT-CASPTRAWATOSROOT,
  "SCRWD-CON-TRAP"=SCRTRADROOSROOT-SCRTRACONOSROOT,
  "SCRWD1R-CON-TRAP"=SCRTRADRORECOSROOT-SCRTRACONOSROOT,
  "SCRWD1R-WD-TRAP"=SCRTRADRORECOSROOT-SCRTRADROOSROOT,
  "SCRWL-CON-TRAP"=SCRTRAWATOSROOT-SCRTRACONOSROOT,
  "SCRWL1R-CON-TRAP"=SCRTRAWATRECOSROOT-SCRTRACONOSROOT,
  "SCRWL1R-WL-TRAP"=SCRTRAWATRECOSROOT-SCRTRAWATOSROOT,
  "SCRPLATE-CON-TRAP"=SCRTRAPLATEOSROOT-SCRTRACONOSROOT,
  levels=design)

# group 2
cont.matrix= makeContrasts(
  "polyAWD-CON"=X35STOTDROOSROOT-X35STOTCONOSROOT,
  "polyAWD1R-CON"=X35STOTDRORECOSROOT-X35STOTCONOSROOT,
  "polyAWD1R-WD"=X35STOTDRORECOSROOT-X35STOTDROOSROOT,
  "polyAWL-CON"=X35STOTWATOSROOT-X35STOTCONOSROOT,
  "polyAWL1R-CON"=X35STOTWATRECOSROOT-X35STOTCONOSROOT,
  "polyAWL1R-WL"=X35STOTWATRECOSROOT-X35STOTWATOSROOT,
  "polyASUB-CON"=X35STOTSUBOSROOT-X35STOTCON2OSROOT,
  "polyASUB1R-CON"=X35STOTSUBRECOSROOT-X35STOTCON2OSROOT,
  "polyASUB1R-SUB"=X35STOTSUBRECOSROOT-X35STOTSUBOSROOT,
  "polyAFIELD-CON"=X35STOTFIELDOSROOT-X35STOTCONOSROOT,
  "polyAFIELD-WAT"=X35STOTFIELDOSROOT-X35STOTWATOSROOT,
  "polyAPLATE-CON"=X35STOTPLATEOSROOT-X35STOTCONOSROOT,
  levels=design)

# group 3
cont.matrix= makeContrasts(
  "RSS1_X35S_CON"=RSS1TRACONOSROOT-X35STRACONOSROOT,
  "RSS1_X35S_CON2"=RSS1TRACON2OSROOT-X35STRACON2OSROOT,
  "QHB_X35S_CON"=QHBTRACONOSROOT-X35STRACONOSROOT,
  "LSI1_X35S_CON"=LSI1TRACONOSROOT-X35STRACONOSROOT,
  "CMZ_X35S_CON"=CMZTRACONOSROOT-X35STRACONOSROOT,
  "CASP_X35S_CON"=CASPTRACONOSROOT-X35STRACONOSROOT,
  "SCR_X35S_CON"=SCRTRACONOSROOT-X35STRACONOSROOT,
  "QHB_RSS1_CON"=QHBTRACONOSROOT-RSS1TRACONOSROOT,
  "LSI1_RSS1_CON"=LSI1TRACONOSROOT-RSS1TRACONOSROOT,
  "CMZ_RSS1_CON"=CMZTRACONOSROOT-RSS1TRACONOSROOT,
  "CASP_RSS1_CON"=CASPTRACONOSROOT-RSS1TRACONOSROOT,
  "SCR_RSS1_CON"=SCRTRACONOSROOT-RSS1TRACONOSROOT,
  "LSI1_QHB_CON"=LSI1TRACONOSROOT-QHBTRACONOSROOT,
  "CMZ_QHB_CON"=CMZTRACONOSROOT-QHBTRACONOSROOT,
  "CASP_QHB_CON"=CASPTRACONOSROOT-QHBTRACONOSROOT,
  "SCR_QHB_CON"=SCRTRACONOSROOT-QHBTRACONOSROOT,
  "CMZ_LSI1_CON"=CMZTRACONOSROOT-LSI1TRACONOSROOT,
  "CASP_LSI1_CON"=CASPTRACONOSROOT-LSI1TRACONOSROOT,
  "SCR_LSI1_CON"=SCRTRACONOSROOT-LSI1TRACONOSROOT,
  "CASP_CMZ_CON"=CASPTRACONOSROOT-CMZTRACONOSROOT,
  "SCR_CMZ_CON"=SCRTRACONOSROOT-CMZTRACONOSROOT,
  "SCR_CASP_CON"=SCRTRACONOSROOT-CASPTRACONOSROOT,
  levels=design)

# group 4
Groups
cont.matrix= makeContrasts(
  "X35SWL-WD-TRAP"=X35STRAWATOSROOT-X35STRADROOSROOT,
  "X35SWL1R-WD-TRAP"=X35STRAWATRECOSROOT-X35STRADROOSROOT,
  "X35SSUB-WD-TRAP"=X35STRASUBOSROOT-X35STRADROOSROOT,
  "X35SSUB1R-WD-TRAP"=X35STRASUBRECOSROOT-X35STRADROOSROOT,
  "X35SFIELD-WD-TRAP"=X35STRAFIELDOSROOT-X35STRADROOSROOT,
  "X35SPLATE-WD-TRAP"=X35STRAPLATEOSROOT-X35STRADROOSROOT,
  "LSI1WL-WD-TRAP"=LSI1TRAWATOSROOT-LSI1TRADROOSROOT,
  "LSI1WL1R-WD-TRAP"=LSI1TRAWATRECOSROOT-LSI1TRADROOSROOT,
  "LSI1FIELD-WD-TRAP"=LSI1TRAFIELDOSROOT-LSI1TRADROOSROOT,
  "RSS1WL-WD-TRAP"=RSS1TRAWATOSROOT-RSS1TRADROOSROOT,
  "RSS1WL1R-WD-TRAP"=RSS1TRAWATRECOSROOT-RSS1TRADROOSROOT,
  "RSS1SUB-WD-TRAP"=RSS1TRASUBOSROOT-RSS1TRADROOSROOT,
  "RSS1SUB1R-WD-TRAP"=RSS1TRASUBRECOSROOT-RSS1TRADROOSROOT,
  "RSS1FIELD-WD-TRAP"=RSS1TRAFIELDOSROOT-RSS1TRADROOSROOT,
  "RSS1PLATE-WD-TRAP"=RSS1TRAPLATEOSROOT-RSS1TRADROOSROOT,
  "QHBFIELD-WD-TRAP"=QHBTRAFIELDOSROOT-QHBTRADROOSROOT,
  "QHBPLATE-WD-TRAP"=QHBTRAPLATEOSROOT-QHBTRADROOSROOT,
  "CMZWL-WD-TRAP"=CMZTRAWATOSROOT-CMZTRADROOSROOT,
  "CMZWL1R-WD-TRAP"=CMZTRAWATRECOSROOT-CMZTRADROOSROOT,
  "CMZPLATE-WD-TRAP"=CMZTRAPLATEOSROOT-CMZTRADROOSROOT,
  "CASPWL-WD-TRAP"=CASPTRAWATOSROOT-CASPTRADROOSROOT,
  "CASPWL1R-WD-TRAP"=CASPTRAWATRECOSROOT-CASPTRADROOSROOT,
  "SCRWL-WD-TRAP"=SCRTRAWATOSROOT-SCRTRADROOSROOT,
  "SCRWL1R-WD-TRAP"=SCRTRAWATRECOSROOT-SCRTRADROOSROOT,
  "SCRPLATE-WD-TRAP"=SCRTRAPLATEOSROOT-SCRTRADROOSROOT,
  levels=design)

#shortName="CELL_TYPE_GH_WL" done
Groups
cont.matrix= makeContrasts(
  "X35SWD1R-WL-TRAP"=X35STRADRORECOSROOT-X35STRAWATOSROOT,
  "X35SSUB-WL-TRAP"=X35STRASUBOSROOT-X35STRAWATOSROOT,
  "X35SSUB1R-WL-TRAP"=X35STRASUBRECOSROOT-X35STRAWATOSROOT,
  "X35SFIELD-WL-TRAP"=X35STRAFIELDOSROOT-X35STRAWATOSROOT,
  "X35SPLATE-WL-TRAP"=X35STRAPLATEOSROOT-X35STRAWATOSROOT,
  "LSI1WD1R-WL-TRAP"=LSI1TRADRORECOSROOT-LSI1TRAWATOSROOT,
  "LSI1FIELD-WL-TRAP"=LSI1TRAFIELDOSROOT-LSI1TRAWATOSROOT,
  "RSS1WD1R-WL-TRAP"=RSS1TRADRORECOSROOT-RSS1TRAWATOSROOT,
  "RSS1SUB-WL-TRAP"=RSS1TRASUBOSROOT-RSS1TRAWATOSROOT,
  "RSS1SUB1R-WL-TRAP"=RSS1TRASUBRECOSROOT-RSS1TRAWATOSROOT,
  "RSS1FIELD-WL-TRAP"=RSS1TRAFIELDOSROOT-RSS1TRAWATOSROOT,
  "RSS1PLATE-WL-TRAP"=RSS1TRAPLATEOSROOT-RSS1TRAWATOSROOT,
  "CMZWD1R-WL-TRAP"=CMZTRADRORECOSROOT-CMZTRAWATOSROOT,
  "CMZPLATE-WL-TRAP"=CMZTRAPLATEOSROOT-CMZTRAWATOSROOT,
  "CASPWD1R-WL-TRAP"=CASPTRADRORECOSROOT-CASPTRAWATOSROOT,
  "SCRWL1R-WL-TRAP"=SCRTRAWATRECOSROOT-SCRTRAWATOSROOT,
  "SCRPLATE-WL-TRAP"=SCRTRAPLATEOSROOT-SCRTRAWATOSROOT,
  levels=design)

#shortName="TREATMENT_GH_SUB" done 
Groups
cont.matrix= makeContrasts(
  "X35SWD1R-SUB-TRAP"=X35STRADRORECOSROOT-X35STRASUBOSROOT,
  "X35SWL1R-SUB-TRAP"=X35STRAWATRECOSROOT-X35STRASUBOSROOT,
  "X35SFIELD-SUB-TRAP"=X35STRAFIELDOSROOT-X35STRASUBOSROOT,
  "X35SPLATE-SUB-TRAP"=X35STRAPLATEOSROOT-X35STRASUBOSROOT,
  "RSS1WD1R-SUB-TRAP"=RSS1TRADRORECOSROOT-RSS1TRASUBOSROOT,
  "RSS1WL1R-SUB-TRAP"=RSS1TRAWATRECOSROOT-RSS1TRASUBOSROOT,
  "RSS1FIELD-SUB-TRAP"=RSS1TRAFIELDOSROOT-RSS1TRASUBOSROOT,
  "RSS1PLATE-SUB-TRAP"=RSS1TRAPLATEOSROOT-RSS1TRASUBOSROOT,
  "X35SWD1R-SUB1R-TRAP"=X35STRADRORECOSROOT-X35STRASUBRECOSROOT,
  "RSS1WD1R-SUB1R-TRAP"=RSS1TRADRORECOSROOT-RSS1TRASUBRECOSROOT,
  "X35SWL1R-SUB1R-TRAP"=X35STRAWATRECOSROOT-X35STRASUBRECOSROOT,
  "RSS1WL1R-SUB1R-TRAP"=RSS1TRAWATRECOSROOT-RSS1TRASUBRECOSROOT,
  "X35SWD1R-WL1R-TRAP"=X35STRADRORECOSROOT-X35STRAWATRECOSROOT,
  "RSS1WD1R-WL1R-TRAP"=RSS1TRADRORECOSROOT-RSS1TRAWATRECOSROOT,
  "CMZWD1R-WL1R-TRAP"=CMZTRADRORECOSROOT-CMZTRAWATRECOSROOT,
  "CASPWD1R-WL1R-TRAP"=CASPTRADRORECOSROOT-CASPTRAWATRECOSROOT,
  "SCRWD1R-WL1R-TRAP"=SCRTRADRORECOSROOT-SCRTRAWATRECOSROOT,
  "LSI1WD1R-WL1R-TRAP"=LSI1TRADRORECOSROOT-LSI1TRAWATRECOSROOT, 
  levels=design)

# group 5
cont.matrix= makeContrasts(
  "RSS1_X35S_FIELD"=RSS1TRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "QHB_X35S_FIELD"=QHBTRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "LSI1_X35S_FIELD"=LSI1TRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "NRAMP_X35S_FIELD"=NRAMPTRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "HMA_X35S_FIELD"=HMATRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "EXP_X35S_FIELD"=EXPTRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "SHR_X35S_FIELD"=SHRTRAFIELDOSROOT-X35STRAFIELDOSROOT,
  "RSS1_SHR_FIELD"=RSS1TRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "QHB_SHR_FIELD"=QHBTRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "LSI1_SHR_FIELD"=LSI1TRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "NRAMP_SHR_FIELD"=NRAMPTRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "HMA_SHR_FIELD"=HMATRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "EXP_SHR_FIELD"=EXPTRAFIELDOSROOT-SHRTRAFIELDOSROOT,
  "QHB_RSS1_FIELD"=QHBTRAFIELDOSROOT-RSS1TRAFIELDOSROOT,
  "LSI1_RSS1_FIELD"=LSI1TRAFIELDOSROOT-RSS1TRAFIELDOSROOT,
  "NRAMP_RSS1_FIELD"=NRAMPTRAFIELDOSROOT-RSS1TRAFIELDOSROOT,
  "HMA_RSS1_FIELD"=HMATRAFIELDOSROOT-RSS1TRAFIELDOSROOT,
  "EXP_RSS1_FIELD"=EXPTRAFIELDOSROOT-RSS1TRAFIELDOSROOT,
  "LSI1_QHB_FIELD"=LSI1TRAFIELDOSROOT-QHBTRAFIELDOSROOT,
  "NRAMP_QHB_FIELD"=NRAMPTRAFIELDOSROOT-QHBTRAFIELDOSROOT,
  "HMA_QHB_FIELD"=HMATRAFIELDOSROOT-QHBTRAFIELDOSROOT,
  "EXP_QHB_FIELD"=EXPTRAFIELDOSROOT-QHBTRAFIELDOSROOT,
  "NRAMP_LSI1_FIELD"=NRAMPTRAFIELDOSROOT-LSI1TRAFIELDOSROOT,
  "HMA_LSI1_FIELD"=HMATRAFIELDOSROOT-LSI1TRAFIELDOSROOT,
  "EXP_LSI1_FIELD"=EXPTRAFIELDOSROOT-LSI1TRAFIELDOSROOT,
  "HMA_NRAMP_FIELD"=HMATRAFIELDOSROOT-NRAMPTRAFIELDOSROOT,
  "EXP_NRAMP_FIELD"=EXPTRAFIELDOSROOT-NRAMPTRAFIELDOSROOT,
  "EXP_HMA_FIELD"=EXPTRAFIELDOSROOT-HMATRAFIELDOSROOT,
  levels=design)


# group 6
cont.matrix= makeContrasts(
  "RSS1_X35S_PLATE"=RSS1TRAPLATEOSROOT-X35STRAPLATEOSROOT,
  "QHB_X35S_PLATE"=QHBTRAPLATEOSROOT-X35STRAPLATEOSROOT,
  "SHR_X35S_PLATE"=SHRTRAPLATEOSROOT-X35STRAPLATEOSROOT,
  "CMZ_X35S_PLATE"=CMZTRAPLATEOSROOT-X35STRAPLATEOSROOT,
  "SCR_X35S_PLATE"=SCRTRAPLATEOSROOT-X35STRAPLATEOSROOT,
  "QHB_RSS1_PLATE"=QHBTRAPLATEOSROOT-RSS1TRAPLATEOSROOT,
  "SHR_RSS1_PLATE"=SHRTRAPLATEOSROOT-RSS1TRAPLATEOSROOT,
  "CMZ_RSS1_PLATE"=CMZTRAPLATEOSROOT-RSS1TRAPLATEOSROOT,
  "SCR_RSS1_PLATE"=SCRTRAPLATEOSROOT-RSS1TRAPLATEOSROOT,
  "SHR_QHB_PLATE"=SHRTRAPLATEOSROOT-QHBTRAPLATEOSROOT,
  "CMZ_QHB_PLATE"=CMZTRAPLATEOSROOT-QHBTRAPLATEOSROOT,
  "SCR_QHB_PLATE"=SCRTRAPLATEOSROOT-QHBTRAPLATEOSROOT,
  "CMZ_SHR_PLATE"=CMZTRAPLATEOSROOT-SHRTRAPLATEOSROOT,
  "SCR_SHR_PLATE"=SCRTRAPLATEOSROOT-SHRTRAPLATEOSROOT,
  "SCR_CMZ_PLATE"=SCRTRAPLATEOSROOT-CMZTRAPLATEOSROOT,
  levels=design)

# group 7
cont.matrix= makeContrasts(
  "RSS1_X35S_DRO"=RSS1TRADROOSROOT-X35STRADROOSROOT,
  "QHB_X35S_DRO"=QHBTRADROOSROOT-X35STRADROOSROOT,
  "LSI1_X35S_DRO"=LSI1TRADROOSROOT-X35STRADROOSROOT,
  "CMZ_X35S_DRO"=CMZTRADROOSROOT-X35STRADROOSROOT,
  "CASP_X35S_DRO"=CASPTRADROOSROOT-X35STRADROOSROOT,
  "SCR_X35S_DRO"=SCRTRADROOSROOT-X35STRADROOSROOT,
  "QHB_RSS1_DRO"=QHBTRADROOSROOT-RSS1TRADROOSROOT,
  "LSI1_RSS1_DRO"=LSI1TRADROOSROOT-RSS1TRADROOSROOT,
  "CMZ_RSS1_DRO"=CMZTRADROOSROOT-RSS1TRADROOSROOT,
  "CASP_RSS1_DRO"=CASPTRADROOSROOT-RSS1TRADROOSROOT,
  "SCR_RSS1_DRO"=SCRTRADROOSROOT-RSS1TRADROOSROOT,
  "LSI1_QHB_DRO"=LSI1TRADROOSROOT-QHBTRADROOSROOT,
  "CMZ_QHB_DRO"=CMZTRADROOSROOT-QHBTRADROOSROOT,
  "CASP_QHB_DRO"=CASPTRADROOSROOT-QHBTRADROOSROOT,
  "SCR_QHB_DRO"=SCRTRADROOSROOT-QHBTRADROOSROOT,
  "CMZ_LSI1_DRO"=CMZTRADROOSROOT-LSI1TRADROOSROOT,
  "CASP_LSI1_DRO"=CASPTRADROOSROOT-LSI1TRADROOSROOT,
  "SCR_LSI1_DRO"=SCRTRADROOSROOT-LSI1TRADROOSROOT,
  "CASP_CMZ_DRO"=CASPTRADROOSROOT-CMZTRADROOSROOT,
  "SCR_CMZ_DRO"=SCRTRADROOSROOT-CMZTRADROOSROOT,
  "SCR_CASP_DRO"=SCRTRADROOSROOT-CASPTRADROOSROOT,
  "RSS1_X35S_DROREC"=RSS1TRADRORECOSROOT-X35STRADRORECOSROOT,
  "LSI1_X35S_DROREC"=LSI1TRADRORECOSROOT-X35STRADRORECOSROOT,
  "CMZ_X35S_DROREC"=CMZTRADRORECOSROOT-X35STRADRORECOSROOT,
  "CASP_X35S_DROREC"=CASPTRADRORECOSROOT-X35STRADRORECOSROOT,
  "SCR_X35S_DROREC"=SCRTRADRORECOSROOT-X35STRADRORECOSROOT,
  "LSI1_RSS1_DROREC"=LSI1TRADRORECOSROOT-RSS1TRADRORECOSROOT,
  "CMZ_RSS1_DROREC"=CMZTRADRORECOSROOT-RSS1TRADRORECOSROOT,
  "CASP_RSS1_DROREC"=CASPTRADRORECOSROOT-RSS1TRADRORECOSROOT,
  "SCR_RSS1_DROREC"=SCRTRADRORECOSROOT-RSS1TRADRORECOSROOT,
  "CMZ_LSI1_DROREC"=CMZTRADRORECOSROOT-LSI1TRADRORECOSROOT,
  "CASP_LSI1_DROREC"=CASPTRADRORECOSROOT-LSI1TRADRORECOSROOT,
  "SCR_LSI1_DROREC"=SCRTRADRORECOSROOT-LSI1TRADRORECOSROOT,
  "CASP_CMZ_DROREC"=CASPTRADRORECOSROOT-CMZTRADRORECOSROOT,
  "SCR_CMZ_DROREC"=SCRTRADRORECOSROOT-CMZTRADRORECOSROOT,
  "SCR_CASP_DROREC"=SCRTRADRORECOSROOT-CASPTRADRORECOSROOT,  
  levels=design)

# group 8
cont.matrix= makeContrasts(
  "RSS1_X35S_WAT"=RSS1TRAWATOSROOT-X35STRAWATOSROOT,
  "LSI1_X35S_WAT"=LSI1TRAWATOSROOT-X35STRAWATOSROOT,
  "CMZ_X35S_WAT"=CMZTRAWATOSROOT-X35STRAWATOSROOT,
  "CASP_X35S_WAT"=CASPTRAWATOSROOT-X35STRAWATOSROOT,
  "SCR_X35S_WAT"=SCRTRAWATOSROOT-X35STRAWATOSROOT,
  "LSI1_RSS1_WAT"=LSI1TRAWATOSROOT-RSS1TRAWATOSROOT,
  "CMZ_RSS1_WAT"=CMZTRAWATOSROOT-RSS1TRAWATOSROOT,
  "CASP_RSS1_WAT"=CASPTRAWATOSROOT-RSS1TRAWATOSROOT,
  "SCR_RSS1_WAT"=SCRTRAWATOSROOT-RSS1TRAWATOSROOT,
  "CMZ_LSI1_WAT"=CMZTRAWATOSROOT-LSI1TRAWATOSROOT,
  "CASP_LSI1_WAT"=CASPTRAWATOSROOT-LSI1TRAWATOSROOT,
  "SCR_LSI1_WAT"=SCRTRAWATOSROOT-LSI1TRAWATOSROOT,
  "CASP_CMZ_WAT"=CASPTRAWATOSROOT-CMZTRAWATOSROOT,
  "SCR_CMZ_WAT"=SCRTRAWATOSROOT-CMZTRAWATOSROOT,
  "SCR_CASP_WAT"=SCRTRAWATOSROOT-CASPTRAWATOSROOT,
  "RSS1_X35S_WATREC"=RSS1TRAWATRECOSROOT-X35STRAWATRECOSROOT,
  "LSI1_X35S_WATREC"=LSI1TRAWATRECOSROOT-X35STRAWATRECOSROOT,
  "CMZ_X35S_WATREC"=CMZTRAWATRECOSROOT-X35STRAWATRECOSROOT,
  "CASP_X35S_WATREC"=CASPTRAWATRECOSROOT-X35STRAWATRECOSROOT,
  "SCR_X35S_WATREC"=SCRTRAWATRECOSROOT-X35STRAWATRECOSROOT,
  "LSI1_RSS1_WATREC"=LSI1TRAWATRECOSROOT-RSS1TRAWATRECOSROOT,
  "CMZ_RSS1_WATREC"=CMZTRAWATRECOSROOT-RSS1TRAWATRECOSROOT,
  "CASP_RSS1_WATREC"=CASPTRAWATRECOSROOT-RSS1TRAWATRECOSROOT,
  "SCR_RSS1_WATREC"=SCRTRAWATRECOSROOT-RSS1TRAWATRECOSROOT,
  "CMZ_LSI1_WATREC"=CMZTRAWATRECOSROOT-LSI1TRAWATRECOSROOT,
  "CASP_LSI1_WATREC"=CASPTRAWATRECOSROOT-LSI1TRAWATRECOSROOT,
  "SCR_LSI1_WATREC"=SCRTRAWATRECOSROOT-LSI1TRAWATRECOSROOT,
  "CASP_CMZ_WATREC"=CASPTRAWATRECOSROOT-CMZTRAWATRECOSROOT,
  "SCR_CMZ_WATREC"=SCRTRAWATRECOSROOT-CMZTRAWATRECOSROOT,
  "SCR_CASP_WATREC"=SCRTRAWATRECOSROOT-CASPTRAWATRECOSROOT,
  levels=design)

# group 9
cont.matrix= makeContrasts(
  "RSS1_X35S_SUB"=RSS1TRASUBOSROOT-X35STRASUBOSROOT,
  "RSS1_X35S_SUBREC"=RSS1TRASUBRECOSROOT-X35STRASUBRECOSROOT,
  levels=design)

# group 10
cont.matrix= makeContrasts(
  "X35S_PLATE_CON"=X35STRAPLATEOSROOT-X35STRACONOSROOT,
  "X35S_FIELD_CON"=X35STRAFIELDOSROOT-X35STRACONOSROOT,
  "X35S_CON2_CON"=X35STRACON2OSROOT-X35STRACONOSROOT,
  "X35S_PLATE_CON2"=X35STRAPLATEOSROOT-X35STRACON2OSROOT,
  "X35S_FIELD_CON2"=X35STRAFIELDOSROOT-X35STRACON2OSROOT,
  "X35S_PLATE_FIELD"=X35STRAPLATEOSROOT-X35STRAFIELDOSROOT,
  "RSS1_PLATE_CON"=RSS1TRAPLATEOSROOT-RSS1TRACONOSROOT,
  "RSS1_FIELD_CON"=RSS1TRAFIELDOSROOT-RSS1TRACONOSROOT,
  "RSS1_CON2_CON"=RSS1TRACON2OSROOT-RSS1TRACONOSROOT,
  "RSS1_PLATE_CON2"=RSS1TRAPLATEOSROOT-RSS1TRACON2OSROOT,
  "RSS1_FIELD_CON2"=RSS1TRAFIELDOSROOT-RSS1TRACON2OSROOT,
  "RSS1_PLATE_FIELD"=RSS1TRAPLATEOSROOT-RSS1TRAFIELDOSROOT,
  "QHB_PLATE_CON"=QHBTRAPLATEOSROOT-QHBTRACONOSROOT,
  "QHB_FIELD_CON"=QHBTRAFIELDOSROOT-QHBTRACONOSROOT,
  "QHB_PLATE_FIELD"=QHBTRAPLATEOSROOT-QHBTRAFIELDOSROOT,
  "CMZ_PLATE_CON"=CMZTRAPLATEOSROOT-CMZTRACONOSROOT,
  "SCR_PLATE_CON"=SCRTRAPLATEOSROOT-SCRTRACONOSROOT,
  "SHR_PLATE_FIELD"=SHRTRAPLATEOSROOT-SHRTRAFIELDOSROOT,
  levels=design)


#### Fit and do Diff Expression
v <- voom(y, design, plot = TRUE,normalize.method ="quantile")
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

# --

### 
dev.off()

## Prepare for gene annotation
annotationAvail <- F
if (annotationAvail){
  cat("Reading annotation file \n")
  genealiasfile <- "gene_aliases.txt"
  ID2Symbol <- getGeneSymbols(genealiasfile)
} else{cat("Annotation file unavailable \n")}


## --

DEList <- list() 

for (contrast in colnames(cont.matrix)){
  print(contrast)
  ## Sorting by none ensures all contrasts will be in the same order
  tmp <- topTable(fit2, coef=contrast,number = Inf,sort.by = "none")
  #
  pValpassed <- table(tmp$adj.P.Val < 0.05)[2]
  cat ("Number of genes with pVal < 0.05 on ",contrast,":",pValpassed,"\n")
  
  
  ## Write genes that are up or downregulated (logFC > 0; logFC < 0)
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


tmpSave <- paste(imgDir,"VolcanoPlots_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")
makeVolcanoPlots(DEList,pValCut=0.01,logCut=2,plotGenes=F) #plotGenes=T to print genes in the plot
dev.off()

### Condense into a single list
## Convert everything to a single table
DE_All <- condenseListTables(DEList) ## Use a custom function

DE_All <- DE_All[,grep("logFC|adj.P.Val",colnames(DE_All),value = T)] #Remove unwanted columns

#### 
DE_All_T <- cbind(DE_All_T,meanTPM[rownames(DE_All),])

## set adj pvalue cutoff
cutoffP=0.01
##function to indicate significant change only (up or down) for each comparison
anno_DEGs=function(df_expr,df_limma,cutoffP){
  ##   
  cutoffP=as.numeric(cutoffP)
  df_expr=as.data.frame(df_expr)
  df_limma=as.data.frame(df_limma)
  adj=grep("adj",colnames(df_limma),value = T)
  fcs=grep("logFC",colnames(df_limma),value = T)
  df_tpm_expr=data.frame(matrix(nrow=nrow(df_expr),ncol = length(fcs)))
  colnames(df_tpm_expr)=paste0(fcs,"_padj<",cutoffP)
  rownames(df_tpm_expr)=rownames(df_expr)
  ## 
  for (i in 1:length(fcs)){
    df_tpm_expr[,i]=""
    for (j in 1:nrow(df_tpm_expr)){
      if ((df_limma[j,fcs[i]]>1) & ((df_limma[j,adj[i]])<cutoffP)){
        df_tpm_expr[j,i]="UP"
      }
      if ((df_limma[j,fcs[i]]<(-1)) & ((df_limma[j,adj[i]])<cutoffP)){
        df_tpm_expr[j,i]="DOWN"
      }
    }
  }
  df_expr=cbind(df_expr,df_tpm_expr)
  return(df_expr)
}
## 
TPM_DEG=anno_DEGs(meanTPM,DE_All,cutoffP = cutoff)
tmpSave <- paste(outDir,"TPM_DEG_Call_",shortName,".csv",sep="")
write.csv(x = TPM_DEG,file = tmpSave,quote = F,row.names = T)
## 
tmpSave <- paste(outDir,"DEG_AllContrasts_",shortName,".csv",sep="")
write.csv(x = DE_All,file = tmpSave,quote = F,row.names = T)

## Close main images
dev.off()
dev.off()
sessionInfo()
sink()

#############################################################################################