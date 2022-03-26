setwd("")
## lists file of meme matrices
All=grep(".meme$",list.files(),value=T)
fimo_function <- function(x){
  setwd("/bigdata/serreslab/mreynoso/Big_Exp/results/DPI/AME/For_FIMO")
  ## lists file of meme matrices
  All=grep(".meme$",list.files(),value=T)
  ## command for FIMO over the whole genome
  file_ds=paste0("module load meme; module load bedtools; fimo --thresh 0.01 --max-stored-scores 1000000000 -o ",Fimo_meme[x],"_FIMO ",Fimo_meme[x],".meme Oryza_sativa.IRGSP-1.0.30.dna.genome.fa")
  ## call command for FIMO
  system(file_ds)
  ## setwd to result folder
  setwd(paste0("/For_FIMO/",Fimo_meme[x],"_FIMO/"))
  ## export bed file of results with locations
  file_ss=paste0("module load bedtools; cut -f3,4,5 fimo.tsv > ", Fimo_meme[x] ,".bed; bedtools sort -i ",Fimo_meme[x],".bed > ",Fimo_meme[x],".sorted.bed; ")
  ## call command for FIMO
  system(file_ss)
  ## libraries
  library(ChIPpeakAnno)
  library(GenomicRanges)
  # reads annotation of genes GRanges
  eByg <- readRDS("Genes.RDS")
  # reads annotation of THSs
  THS_gr <- readRDS("ATAC_GH_FIELD_35S.RDS")
  # reads differential THSs
  Diff_THS_gr=readRDS("Diff_THS_gr.RDS")
  # reads annotation of genes
  annotation=readRDS("annotation_040620.RDS")
  setwd(paste0("/For_FIMO/",Fimo_meme[x],"_FIMO/"))
  ## reads FIMO results  
  file_fimo=read.delim(paste0("fimo.tsv"),comment.char = "#")
  ## adds colnames
    colnames(file_fimo)=gsub("sequence_name","seqnames",colnames(file_fimo))
    ## converts to Granges
    file_fimo_gr=GRanges(file_fimo)
    ## annotate with respect to genes
    annoGR=as.data.frame(annotatePeakInBatch(file_fimo_gr,eByg))
    ## find overlaps with THSs locations
    object=findOverlaps(file_fimo_gr,THS_gr)
    ## aux dataframe with overlaps
    df=as.data.frame(object)
    ## Add names to target regions
    names(file_fimo_gr)=1:length(file_fimo_gr)
    ## find overlaps with dTHSs locations
    object2=findOverlaps(file_fimo_gr,Diff_THS_gr)
    ## aux dataframe with overlaps
    dfDEG=as.data.frame(object2)
    ## addition to annotated THSs
    annoGR$THS=df$subjectHits[match(annoGR$peak,df$queryHits)]
    ## addition of labels
    annoGR$labels=Diff_THS_gr$labels[dfDEG$subjectHits[match(annoGR$peak,dfDEG$queryHits)]]
    ## addition to dTHSs
    annoGR$dTHS=dfDEG$subjectHits[match(annoGR$peak,dfDEG$queryHits)]
    ## addition of description
    annoGR$description=annotation$Description[match(annoGR$feature,annotation$GENEID)]
    ## addition of short names
    annoGR$CGSNL.Gene.Symbol=annotation$CGSNL.Gene.Symbol[match(annoGR$feature,annotation$GENEID)]
    ## filter step to upstream closer than 2kb, either inside or downstream less than 1kb
    filtered_annoGR=annoGR[(annoGR$insideFeature=="upstream"&annoGR$shortestDistance<2000)|(annoGR$insideFeature=="downstream"&annoGR$shortestDistance<1000)|(annoGR$insideFeature=="inside"),]
    write.csv(filtered_annoGR,paste0("annoGR_",Fimo_meme[x],".csv"))
}
#### Run on cluster or individually locally 
library(BiocParallel); library(BatchJobs)
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="4:00:00", ntasks=length(Fimo_meme), ncpus=1, memory="8G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(Fimo_meme), resources=resources, cluster.functions=cluster.functions)
register(param)
Funct <- bplapply(seq(along=Fimo_meme), fimo_function) # Check status with qstat
