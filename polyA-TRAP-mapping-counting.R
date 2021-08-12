##################################################################################
############ Mapping polyA and TRAP RNA-seq reads to rice genome #################
##################################################################################

##### Fastq files as obtained from the sequencing facility ###
## Loading of libraries required throughout the procedure ###
library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
## Set working directory
setwd("MyDirectory/")
## File containing the location of fastq.gz files and a column for names
targets=read.delim("targets.txt", header=T,comment.char = "#")
input=targets$FileName
library(BiocParallel); library(BatchJobs)
#function definition
fmap <- function(x) {
  ## Several packages needed 
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
  library(BiocParallel); library(BatchJobs)
  setwd("MyDirectory/")
  targets=read.delim("targets.txt", header=T,comment.char = "#")
  input=targets$FileName
  name=targets$SampleName
  ## reference
  reference <- "/Mydirectory/Oryza_sativa.IRGSP-1.0.30.dna.genome"
  ## gff3 genome annotation file
  gff3<-"/Mydirectory/Oryza_sativa.IRGSP-1.0.30.gff3"
  ## command for mapping to the genome
  tophat_command <- paste0("module load bowtie2; module load tophat; tophat -p 8 -g 1 -i 50 -I 3000 -G ",gff3," -o ", name[x], ".tophat2 ", reference, " ", input[x])
  ## execution in the system
  system(tophat_command)
  ## sort bam files in R
  sortBam(file=paste(name[x], ".tophat2/accepted_hits.bam", sep=""), destination=paste(name[x], ".tophat2/accepted_hits", sep="")) 
  ## index bam files in R
  indexBam(paste(name[x], ".tophat2/accepted_hits.bam", sep=""))
}
## Call to the function to execute in multiple cores
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(input), ncpus=length(input), memory="20G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input), resources=resources, cluster.functions=cluster.functions)
register(param)
bplapply(seq(along=input), fmap)
## Function fmap can be called with other loop to obtain mapped reads
##########################################################

##########################################################
########## Read counts over exonic regions ###############
##########################################################

## Load R
R
## Set working directory
setwd("/Mydirectory/results")
## File containing the location of fastq.gz files and a column for names
targets=read.delim("/Mydirectory/targets.txt", header=T,comment.char = "#")
## Naming code = Promoter + Sample type + Species + Tissue + Replicate number

fcount <- function(x) {
  ## Loading of libraries
  library(systemPipeR); library(BiocParallel); library(GenomicAlignments);library(GenomicFeatures);library(systemPipeR)
  setwd("/Mydirectory/results")
  targets=read.delim("/Mydirectory/targets.txt", header=T,comment.char = "#")
  ## Name of the samples
  name=targets$SampleName
  ## location of bam files
  mybams <- paste(name, ".tophat2/accepted_hits.bam", sep="")                                                                                                                                                                    
  ## Name assigned to the bam vector
  names(mybams)=name
  ## Create object with list of bam files
  bfl <- BamFileList(as.character(mybams), yieldSize=50000, index=character())
  ## Assigned to the bam list
  names(bfl)=names(mybams)
  ## Load of transcript database generated from gff3 using systempipeR function makeTxDbFromGFF over Oryza_sativa.IRGSP-1.0.30.gff3
  txdb <- loadDb("/Mydirectory/txdb.sqlite") 
  ## Functions to extract exons by genes object
  eByg <- exonsBy(txdb, by="gene")
  ## Function to count each bam file over features. Notice data is not strand specific 
  summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=T, inter.feature=FALSE, singleEnd=TRUE)
}
## Call to the function to execute in multiple cores
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(input), ncpus=length(input), memory="1G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input), resources=resources, cluster.functions=cluster.functions)
register(param)
counteByg <- bplapply(seq(along=input), fcount) # Check status with qstat
## Function fcount can be called with other loop to obtain mapped reads

## Assign counts to create a new dataframe 
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## Assign rownames corresponding to each gene and colnames corresponding to each sample 
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- name
## Export results
write.table(countDFeByg, "/Mydirectory/countExonsBygenes.xls", quote=FALSE, sep="\t", col.names = NA)
##########################################################

