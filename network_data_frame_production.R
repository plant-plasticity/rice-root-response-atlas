### Produce dataframe for networks on Cyverse ### 
##working directory
setwd("")
### Group of genes of interest
Patterns=read.csv("Patterns.csv")
### AME Result filtered by expression level
AME_THS_ann_filt_level=read.csv("AME_filtered.csv")

### Motifs enriched
enmotif=c(unique(as.character(AME_THS_ann_filt_level$CIS.BP.DB.Motif.ID)))
### avoid NA
enmotif=enmotif[!is.na(enmotif)]

### Gene targets
GeneTargets=Patterns$GENEID[Patterns$Pattern%in%"Pattern_4"]
### avoid NA
GeneTargets=GeneTargets[!is.na(GeneTargets)]
### First add motif 1 
  test=file.exists(paste0("For_FIMO/",enmotif,".meme_FIMO/annoGR_",enmotif,".csv"))|file.exists(paste0("For_FIMO/",enmotif,".meme_FIMO/annoGR_",enmotif,".meme.csv"))
### Data frame with information about the motif
  df=data.frame(enmotif,test,AME_THS_ann_filt_level$TF.Fam[match(enmotif,AME_THS_ann_filt_level$CIS.BP.DB.Motif.ID)])
  colnames(df)[3]="TF.Fam"
### Motif similarities bases
  dfLessStrict=readRDS("dfLessStrict.RDS")
  dfMiddleStrict=readRDS("dfMiddleStrict.RDS")
  dfMoreStrict=readRDS("dfMoreStrict.RDS")
### Motif similarities bases  
  df$lessStrict=0
  df$middleStrict=0
  df$moreStrict=0
  for (i in 1:length(df$enmotif)){
    ## rest of motifs
    df$lessStrict[i]=grep(df$enmotif[i],dfLessStrict$Group)[1]
    df$middleStrict[i]=grep(df$enmotif[i],dfMiddleStrict$Group)[1]
    df$moreStrict[i]=grep(df$enmotif[i],dfMoreStrict$Group)[1]
  }
  ### 
  write.csv(df,paste0("Stats_for_motifs.csv"))
  n=match("TRUE",test)
  ## load fimo annotated FIMO results for first enriched motif
  if (file.exists(paste0("For_FIMO/",enmotif[n],".meme_FIMO/annoGR_",enmotif[n],".csv"))){
    mot_targets=read.csv(paste0("For_FIMO/",enmotif[n],".meme_FIMO/annoGR_",enmotif[n],".csv"))
    mot_targets=mot_targets[mot_targets$feature%in%GeneTargets,1:ncol]
  }
  if ((n+1)<length(enmotif)){
    ## load fimo annotated FIMO results for each enriched motif
    for (s in (n+1):length(enmotif)){
      if (file.exists(paste0("For_FIMO/",enmotif[s],".meme_FIMO/annoGR_",enmotif[s],".csv"))){
        mot_targets2=read.csv(paste0("For_FIMO/",enmotif[s],".meme_FIMO/annoGR_",enmotif[s],".csv"))  
        mot_targets2=mot_targets2[mot_targets2$feature%in%GeneTargets,1:ncol]
        mot_targets=rbind(mot_targets,mot_targets2)
      }
    
  }
  }
## Filter location
  mot_targets=mot_targets[!((mot_targets$insideFeature=="upstream"&mot_targets$shortestDistance>2000|mot_targets$insideFeature=="downstream"&mot_targets$shortestDistance>1000)),]
## Filter inside THS only
  mot_targets=mot_targets[!mot_targets$type=="p",]
## add gene ID
  mot_targets=cbind(mot_targets,Patterns[match(mot_targets$feature,Patterns$GENEID),])
## 
  mot_targets$TF.Fam=df$TF.Fam[match(mot_targets$motif_alt_id,df$enmotif)]

# Export network
write.csv(mot_targets[,c(8,14,22,24:27,29,36:39,41,105:121)],paste0("_Network_","Pattern_4","_simplified_220125.csv"))   
if(length(mot_targets[,1])>0){
  ### collapsing similar motifs
  library(GenomicRanges)
  ## targets each group collapsing based on overlap of locations ##
  mot_targets$lessStrict=df$lessStrict[match(mot_targets$motif_alt_id,df$enmotif)]
  mot_targets$middleStrict=df$middleStrict[match(mot_targets$motif_alt_id,df$enmotif)]
  mot_targets$moreStrict=df$moreStrict[match(mot_targets$motif_alt_id,df$enmotif)]
  ## make GRanges
  mot_targetsgr=GRanges(mot_targets)
  ## reduce Granges
  mot_targetsgrRed=reduce(mot_targetsgr)
  ## find overlaps over GRanges 
  rep=data.frame(findOverlaps(mot_targetsgr,mot_targetsgrRed))
  ## Use overlaps to find collapsing 
  mot_targets$Collapse=rep$subjectHits[match(1:length(mot_targetsgr),rep$queryHits)]
  mot_targets$Coll_Strict=paste0(mot_targets$moreStrict,"_",mot_targets$Collapse)
  mot_targets$Coll_Less_Strict=paste0(mot_targets$lessStrict,"_",mot_targets$Collapse)
  ## 
  collapse_list=unique(mot_targets$Coll_Strict)
  ## export collapsed for strict or less similarity
  collapse_list_less=unique(mot_targets$Coll_Less_Strict)
  mot_targets_collapsed_strict=mot_targets[match(collapse_list,mot_targets$Coll_Strict),]
  mot_targets_collapsed_less=mot_targets[match(collapse_list_less,mot_targets$Coll_Less_Strict),]
## export collapsed for strict or less similarity
  write.csv(mot_targets_collapsed_strict[,],paste0("Patterns_Network_","_THSs_collapsed_strict.csv"))
  write.csv(mot_targets_collapsed_less[,],paste0("Patterns_Network_","_THSs_collapsed_loose.csv"))
}
