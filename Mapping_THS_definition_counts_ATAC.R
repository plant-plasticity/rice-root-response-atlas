##################################################################################
############ Mapping ATAC-seq reads to rice genome #################
##################################################################################
## Load R
R
## Set working directory
setwd("/mydirectory/")
## Loading of libraries
library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
library(BiocParallel); library(BatchJobs)
## File containing the location of paired fastq.gz files and a column for names
targets=read.delim("/mydirectory/targetsATAC.txt", header=T,comment.char = "#")
## vector with location for pair 1
input1=targets$FileName1
## vector with location for pair 2
input2=targets$FileName2
#function definition
fmap <- function(x) {
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
  library(BiocParallel); library(BatchJobs)
  setwd("/mydirectory/")
  ## File containing the location of paired fastq.gz files and a column for names
  targets=read.delim("/mydirectory/targetsATAC.txt", header=T,comment.char = "#")
  ## vector with location for pair 1
  input1=targets$FileName1
  ## vector with location for pair 2
  input2=targets$FileName2
  ## vector with sample name
  name=targets$SampleName
  ## vector with output bam name
  bams<-paste0(targets[,3],".filt.bam")
  ## vector with output bai name
  bais<-paste0(targets[,3],".filt")
  ## Indexed name for chloroplast genome
    genomechl<-"OsChloroplast"
    ## Location for chloroplast genome
    referencechl <- paste0("/mydirectory/",genomechl)
    ## vector aligned to chloroplast genome
    bams_chl=paste0(name,"_",genomechl,".bam")
    ## Filter command to keep unaligned reads to chloroplast
    bowtie_filter_command=paste0("module load samtools; module load bowtie2; bowtie2 --local -p 8 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x ",referencechl," --un-conc ",name[x], "_" ,genomechl, ".unaligned -1 ", input1[x]," -2 ",input2[x]," | samtools view -bS - > ",bams_chl[x])
    ## Filter command 
    system(bowtie_filter_command)
    ## Vector of filtered reads pair 1
    unalign1chl=paste0(name, "_" ,genomechl, ".1.unaligned")
    ## Vector of filtered reads pair 2
    unalign2chl=paste0(name, "_" ,genomechl, ".2.unaligned")
    ## Indexed name for mitochondrial genome
    genomemit<-"OsMitochondria"
    ## Location for mitochondrial genome
    referencemit <- paste0("/mydirectory/",genomemit)
    ## vector aligned to mitochondrial genome
    bamsmit=paste0(name,"_",genomemit,".bam")
    ## Filter command keep unaligned reads reads to chloroplast and mitochondria
    bowtie_mit_command=paste0("module load samtools; module load bowtie2; bowtie2 --local -p 8 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x ",referencemit," --un-conc ",name[x], "_" ,genomemit, ".unaligned -1 ",unalign1chl[x]," -2 ",unalign2chl[x], " | samtools view -bS - > ",bamsmit[x])
    ## Filter command 
    system(bowtie_mit_command)
    ## Vector of filtered reads to chloroplast and mitochondria pair 1
    unalign1chlmit=paste0(name, "_" ,genomemit, ".1.unaligned")
    ## Vector of filtered reads to chloroplast and mitochondria pair 1
    unalign2chlmit=paste0(name, "_" ,genomemit, ".2.unaligned")
    ## Reference genome
    reference <- "/mydirectory/Oryza_sativa.IRGSP-1.0.30.dna.genome"
    ## Reference annotation
    gff3<-"/mydirectory/Oryza_sativa.IRGSP-1.0.30.gff3"
    ## Mapping command to keep alignments with Q>10
    mapping_command <- paste0("module load samtools; module load bowtie2; bowtie2 -p 8 -x ",reference," -1 ", unalign1[x]," -2 ",unalign2[x]," | samtools view -bS -q 10  - > ",bams[x])
    system(mapping_command)
    ## Sort and index results
    sortBam(file=bams[x], destination=bais[x]) 
    indexBam(bams[x])
  
}
## Call to the function to execute in multiple cores
cluster.functions <- makeClusterFunctionsSLURM("mydirectory/slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(input1), ncpus=length(input1), memory="10G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input1), resources=resources, cluster.functions=cluster.functions)
register(param)
bplapply(seq(along=input1), fmap)
## Function fmap can be called with other loop to obtain mapped reads
##################################################################################


##################################################################################
############ THSs definition HOMER #################
##################################################################################
## Set working directory
setwd("/mydirectory/")
## File containing the location of paired fastq.gz files and a column for names
targets=read.delim("targetsATAC.txt", header=T,comment.char = "#")
## vector with sample name
name=targets$SampleName
## bam files vector
mybams <- paste0("/mydirectory/",targets[,3],".filt.bam")
## bam files list
bfl <- BamFileList(as.character(mybams), yieldSize=50000, index=character())                                                                                                                  
names(bfl)=name
## names of THS final file
mer=paste0(targets$SampleName,".q10.homer.region.mD150.file.merged.bed")
## function definition
fhomer=function(x){
  ## packages required
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead)
  library(BiocParallel); library(BatchJobs) 
  setwd("/mydirectory/")
  ## File containing the location of paired fastq.gz files and a column for names
  targets=read.delim("targetsATAC.txt", header=T,comment.char = "#")
  ## vector with sample name
  name=targets$SampleName
  ## bam files vector
  mybams <- paste0("/mydirectory/",targets[,3],".filt.bam")
  ## bam files list
  bfl <- BamFileList(as.character(mybams), yieldSize=50000, index=character())                                                                                                                  
  names(bfl)=name
  ## names of THS final file
  mer=paste0(targets$SampleName,".q10.homer.region.mD150.file.merged.bed")
  ## homer command parameters -gsize 3.7e8  -minDist 150
  homer_command=paste0("module load homer; module load samtools; module load bedtools; makeTagDirectory ", names(bfl)[x],".q10.homer.tags/ ", mybams[x],"; findPeaks ",names(bfl)[x],".q10.homer.tags/ -o ",names(bfl)[x],".q10.homer.region.mD150 -gsize 3.7e8  -minDist 150 -region; pos2bed.pl ",names(bfl)[x],".q10.homer.region.mD150 | bedtools sort|bedtools merge > ",names(bfl)[x],".q10.homer.region.mD150.file.merged.bed")
  ## homer command callout
  system(homer_command)
}
## Call to the function to execute in multiple cores
funs <- makeClusterFunctionsSLURM("/bigdata/serreslab/mreynoso/Big_Exp/slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(mybams), ncpus=1, memory="10G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(mybams), resources=resources, cluster.functions=funs)
bplapply(seq(along=mybams), fhomer)
## Function fmap can be called with other loop to obtain THS for each sample
##################################################################################

##################################################################################
############ THSs in more than one replicate and annotatation #################
##################################################################################
## load libraries
library(ChIPpeakAnno);library(GenomicFeatures);library(GenomicRanges)
## File containing the location of paired fastq.gz files and a column for names
targets=read.delim("targetsATAC.txt",header=T)
## Select samples  X35SATCON
targets_sel=targets[targets[,3]%in%grep("X35SATCON",as.character(targets[,3]),value=TRUE),]
## Creates an empty list with the genomic ranges of each THS 
gr35S=GRangesList()
## Loading each set of peaks to the list
for(x in seq(along=targets_sel[,3])) {
  ## Loading as a data frame bed file from homer
  Peak_set<-read.delim(paste0( "./mydirectory/",as.character(targets_sel[x,3]),".q10.homer.region.mD150.file.merged.bed"), header=F)
  # assigning names to columns for granges
  colnames(Peak_set)=c("seqnames","start","end")
  # Converting THSs to granges objects and assigning to the granges list
  gr35S[[x]] <- makeGRangesFromDataFrame(Peak_set, keep.extra.columns=TRUE)
  # Add replicate names to the granges list
  names(gr35S)[x]=paste0(as.character(targets_sel[x,3]))
}
# Validate intersection with overlap of at least 50 bp and merging peaks that overlap 
overlap_35S_merge=findOverlapsOfPeaks(gr35S,minoverlap = 50,connectedPeaks = "merge")
## Save overlap objects 
saveRDS(overlap_35S_merge,"./Lists_overlap_merge_35S_SB_CON_homer.RDS")
## GRanges object for merged peaks of 35S in 2 or more replicates
merged_35S_CON=overlap_35S_merge$mergedPeaks
## Merged THSs for different conditions
gr35S=c(merged_35S_CON,merged_35S_CON2,merged_35S_FIELD,merged_35S_DRO,merged_35S_DRO_REC,merged_35S_SUB,merged_35S_SUB_REC,merged_35S_WAT,merged_35S_WAT_REC)
## Save RDS
saveRDS(gr35S,"./35S_THSs.RDS")
## annotate peaks using txdb transcript database 
annotatedPeak <- annotatePeakInBatch(gr35SFIred, AnnotationData=genes(txdb))
## Create data.frame
df_anno_GH_gr35S <- data.frame(as.data.frame(annotatedPeak))
## export data.frame 
write.csv(df_anno_GH_gr35S,"./df_anno_GH_FI.csv")
##################################################################################

##################################################################################
############ Count ATAC reads over region of interest #################
##################################################################################
library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
library(BatchJobs);library(stringr)
##
setwd("/mydirectory/")
## 
fil=list.files()
##
bams=grep("filt.bam$",fil,value=T)
##
mybams=grep("X35SAT",bams,value=T)
##
bfl <- BamFileList(as.character(mybams), yieldSize=50000, index=character())                                                                                                                   
##
for (i in 1:length(mybams)){
  ##
 names(bfl)[i]=substr(mybams[i],1,(str_locate(mybams[i],".filt.bam$")[1]-1))
}
##
THS<-readRDS("35S_THSs.RDS")
##
counteByg=summarizeOverlaps(THS, bfl[1:length(bfl)], mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=FALSE)
##
cA=assays(counteByg[1:length(THS)])$counts
##
write.csv(cA,"THS_counts.csv")
##################################################################################

##################################################################################
############ Bed files to regions 2kb upstream ATG #################
##################################################################################
## count over 1kb promoter, or more 2kb!? from TSS or ATG?
# generate counts for both 

## load libraries
library(GenomicFeatures)
## load trasncript database
txdb=loadDb("/mydirectory/txdb.sqlite")
## extract CDS regions
cds_os <- cdsBy(txdb, by="tx", use.names=TRUE)
## change names to genes instead of transcripts for rice RAP ID G instead of T
names(cds)=gsub("T","G",names(cds))
## for keeping all genes
cds_subset=cds_os
## Name ouptup name
name_output="Promoter_ATG"
## Indicate the folder containing the file  ##########
wd="/mydirectory/"
####### Change directory #############
setwd(wd)
## select the distance to ATG to take as promoter
disATG=2000
## Creates an object that will contain the bed information (6 columns)
prom_os=as.data.frame(matrix(nrow=length(names(cds_subset)),ncol=6))
## add names to each columns
colnames(prom_os)=c("seqnames","start","end","Tx","score","strand")
prom_os$score="."
## Loop to fill in the information of coordinates
for (j in 1:length(names(cds_subset))){
  ## starts
  x=as.data.frame(start(cds_subset[names(cds_subset)==names(cds_subset)[j]]))
  ## ends
  y=as.data.frame(end(cds_subset[names(cds_subset)==names(cds_subset)[j]]))
  ## strands
  z=as.data.frame(strand(cds_subset[names(cds_subset)==names(cds_subset)[j]]))
  ## names
  prom_os$Tx[j]=z$group_name[1]
  ## chromosome
  prom_os$seqnames[j]=as.character(as.data.frame(seqnames(cds_subset[names(cds_subset)==names(cds_subset)[j]]))$value[1])
  prom_os$strand[j]=as.character(z$value[1])
  if (z$value[1]=="+"){
    prom_os$start[j]=x$value[1]-disATG
    prom_os$end[j]=x$value[1]
  }
  if (z$value[1]=="-"){
    prom_os$start[j]=y$value[1]
    prom_os$end[j]=y$value[1]+disATG
  }
}

## Excluding those that are less than 0 coordinates.
prom_os=prom_os[prom_os$start>0,]
## Option instead of excluding incorporate a 0 instead :
prom_os$start[prom_os$start<0]=0
## Write a bed file
write.table(prom_os,paste0(name_output,disATG,".bed"),quote=F, sep="\t", row.names=F, col.names=F)
################# To select group of genes #######################
#### read the file containing genes of interest
up=read.csv(name_file)
## Indicate which Gene ID are found in the list 
select=match(substr(unique(up$GENE_ID),1,13),substr(names(cds_os),1,13))
## A keep option
select=select[!is.na(select)]
## From the total cordinates of CDS subset it to genes of interest 
cds_subset=cds_os[select]
## Run loop and export
##################################################################################