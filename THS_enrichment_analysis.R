################### Test enrichment #####################
#### read group of genes to test for enrichment 
GeneSet=read.csv("GeneSet.csv")
### Group of unique genes
names_groups=as.character(unique(GeneSet$Pattern))
###  Saves object 
saveRDS(names_groups,"GeneSet.RDS")

##### File of genes and associated THSs ######
dTHSRNA=read.csv("THS_FC_annotated.csv")

### Fasta sequences of THSs #####
library(seqinr)
## fasta 
OsTHS=read.fasta("os_THS.fa")

### Function to select which THSs correspond are nearby genes of interest to export ######
SelectGenesDPI=function(DPI_table,fasta,DPI,THSRNA){
  ## Genes in Geneset
  Enrich=DPI_table$Gene[DPI_table$Pattern==DPI]
  ## THSs ID associated to genes in GeneSet
  Enrich_THS=THSRNA$THS_ID[THSRNA$feature%in%Enrich]
  ## Select which THSs to keep
  selEnr=(names(fasta)%in%Enrich_THS)
  ## return object specifying THS of interest
  return(selEnr)
}
## test of function
selection=SelectGenesDPI(GeneSet,OsTHS,names_groups[[1]],dTHSRNA)
## 
table(selection)
### Function to export fasta files ####
ExportSetFasta=function(selection,fasta,name){
  ## Export fasta derived from fasta files
  write.fasta(fasta[selection],names(fasta)[selection],file.out =paste0(name,".fa"))
}
## setwd directory ###
setwd("")
## application of the functions for all the GeneSets ##
for (i in 1:length(names_groups)){
  ## application of the functions for all the GeneSets ##
  selEnr=SelectGenesDPI(GeneSet,OsTHS,names_groups[[i]],dTHSRNA)
  ## control for each set of genes that contains the rest of THSs universe under this condition ###
  selsetCon=!selEnr
  ##  Export fasta files containing sequences THSs for genes nearby
  ExportSetFasta(selEnr,OsTHS,paste0("Enriched.GeneSet.THS.",names_groups[[i]]))
  ##  Export fasta files containing control sequences THSs 
  ExportSetFasta(selsetCon,OsTHS,paste0("Control.GeneSet.THS.",names_groups[[i]]))
}

############################################
#### AME test comparing Enriched set vs Control set
###############################
## 

setwd("")

names_groups=readRDS("GeneSet.RDS")
## Function to calculate enrichment in sequences set compared to control
fenrich <- function(x){
  ## directory
  setwd("")
  names_groups=readRDS("GeneSet.RDS")
  ###output
  dir=paste0("ame_Enr_GeneSet_20_THS_",names_groups[[x]])
  ## set of THSs of interest
  fas=paste0("Enriched.GeneSet.THS.",names_groups[[x]],".fa")
  ## set of control THSs
  cont=paste0("Control.GeneSet.THS.",names_groups[[x]],".fa")
  ## run command for AME enrichment on CIS BP rice transcription factors and Arabidopsis
  Ame=paste0("module load meme/5.0.2; ame -oc ",dir," --control ",cont,"  --hit-lo-fraction 0.01 --evalue-report-threshold 10 ",fas," OryzaSativaCisBPDBMotifsPlusHRPE.meme ArabidopsisDAPv1.meme")
  system(Ame)
}
#### Cluster running for parallel running. Not exclusive function can be run once at a time.
library(BiocParallel); library(BatchJobs)
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="10:00:00", ntasks=length(names_groups), ncpus=1, memory="12G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(names_groups), resources=resources, cluster.functions=cluster.functions)
register(param)
Function <- bplapply(seq(along=names_groups), fenrich) # Check status with qstat
####################

####### select and save table #########
## collapsing results for different sets
setwd("")
names_groups=readRDS("GeneSet_20_DPI.RDS")
x=2
##
dir=paste0("ame_Enr_GeneSet_20_THS_",names_groups[[x]])
## read results 
AME=read.delim(paste0(dir,"/ame.tsv"))
## Add group names to the result table
AME$Group_name=names_groups[[x]]
## cycle for the rest of genesets
for (x in 1:length(names_groups)){
  dir=paste0("ame_Enr_GeneSet_20_THS_",names_groups[[x]])
  AME2=read.delim(paste0(dir,"/ame.tsv"))
  if(length(AME2[1,])>1){
    AME2$Group_name=names_groups[[x]]
    AME=rbind(AME,AME2)  
  }
}
## export data
write.csv(AME[!is.na(AME$adj_p.value),],paste0("AME_results_THS_","GeneSet_20.csv"))

####### heatmap to visualize enrichment #########
RNATPM=read.csv("TPM_DEG_Call_mar_root_GH_Field_Plate_polyA_TRAP.csv")

## function to add annotation
addanno=function(AME_SUB_pearson,set){
  ## convert gene IDs
  AME_SUB_pearson$GeneID=gsub("g","G",gsub("s","S",substr(AME_SUB_pearson$motif_ID,1,12)))
  ## associate genes with TPM data
  AME_SUB_pearson[,21]=RNATPM[match(AME_SUB_pearson$GeneID,RNATPM$X),set]
  ## read annotation
  annoTF=read.csv("annotation.csv")
  ## add annotation matching gene ID
  AME_SUB_pearson[,22:28]=annoTF[match(AME_SUB_pearson$GeneID,annoTF$Gene),2:8]
  ## return data
  return(AME_SUB_pearson)
}

## read AME results
AMEGroup=read.csv("AME_results_THS_GeneSet_20.csv")
## group of samples to add mean TPM values for 
set="X35STOTCONOSSHOOTB"
## group of 
AMEGroup_ann=addanno(AMEGroup,set)
## group of 
AMEGroup_ann$Class=paste0(AMEGroup_ann$CIS.BP.DB.Motif.ID,AMEGroup_ann$Group_name)

## CisBP comparison of transcription factor domains in rice CisBP
dfLessStrict=readRDS("dfLessStrict.RDS")
dfMoreStrict=readRDS("dfMoreStrict.RDS")

##
unique(AMEGroup_ann$Class)
####### Keep only one group of TFs for group of genes to avoid repetitions
AMEGroup_ann_simpl=AMEGroup_ann[match(unique(AMEGroup_ann$Class),AMEGroup_ann$Class),]
## Add groups for less strict
AMEGroup_ann_simpl$lessStrict=0
AMEGroup_ann_simpl$moreStrict=0

for (i in 1:length(AMEGroup_ann_simpl$moreStrict)){
  AMEGroup_ann_simpl$lessStrict[i]=grep(AMEGroup_ann_simpl$CIS.BP.DB.Motif.ID[i],dfLessStrict$Group)[1]
  AMEGroup_ann_simpl$moreStrict[i]=grep(AMEGroup_ann_simpl$CIS.BP.DB.Motif.ID[i],dfMoreStrict$Group)[1]
}
## export file
write.csv(AMEGroup_ann,"AME_results_THS_GeneSet_20_annotated.csv")

## add library documents
library(ggplot2)
library(viridis)
## Distance for heatmap
X=0.2
class(X)="unit"
## ggplot call for heatmap
pdf("Heatmap_TF_enrichment_THS_GeneSet_20_THS_color_simpl_Spacing_02.pdf")
(ggplot(AMEGroup_ann_simpl,aes(x=Fig2_Pattern,y=GeneID,fill=-log10(adj_p.value))) + 
    theme(strip.text.x = element_text(size = 8, colour = "orange", angle = 90))+
    geom_tile() +
    theme(panel.spacing = X)+
    scale_fill_gradientn(colours=viridis(10,direction=-1),values = c(0.0005,0.5,1),na.value = "black")+
    facet_grid(AMEGroup_ann_simpl$TF.Fam ~ .,scales="free_y",space="free",switch = "y")+
    theme(strip.text.y.left = element_text(angle = 0,size=8),strip.text.x = element_text(angle = 45))
  + theme(axis.text.x = element_text(angle = 45,size = 8,hjust = 1))+
    theme(axis.text.y = element_text(size = 2)))
dev.off()
