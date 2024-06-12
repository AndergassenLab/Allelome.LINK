##*************************************************##
##  Allelome.PRO / Link generate coding info file  ## 
##*************************************************##
## Tim Hasenbein
## Last modification 12.2022
## Creation: 12.2022
## Generate coding info file from RefSeq for Allelome.PRO / Link in BED6 format
## 1. Chromosome, 2. Name, 3. Strand, 4. coding_noncoding
library(plyr)


######------ set environment -----######
input <- "/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/Y_references/annotation/20180203_RefSeq/"
output <- "/Users/timhasenbein/Desktop/"
chr <- c(paste(rep("chr",19),1:19,sep=""),"chrX","chrY")


######------ get reference data -----######
Refseq_IDs <- read.table(file=paste0(input,"20180203_RefSeq_IDs.txt"),sep="\t",header=F,stringsAsFactors=FALSE)
colnames(Refseq_IDs) <- c("transcript_id","chr","strand","start","end","name")
Refseq_IDs <- Refseq_IDs[Refseq_IDs$chr %in% chr,]
# extract NM, XM = coding and NR and XR = noncoding
Refseq_IDs$coding_info <- apply(Refseq_IDs,1,function(x){
  return(unlist(strsplit(x["transcript_id"],split = "_"))[1])
})
Refseq_IDs[Refseq_IDs$coding_info %in% c("NM","XM"),"coding_noncoding"] <- "coding"
Refseq_IDs[Refseq_IDs$coding_info %in% c("NR","XR"),"coding_noncoding"] <- "noncoding"
# if any isoform is coding, coding is assigned. if all noncoding, noncoding is assigned
temp <- NULL
temp <- ddply(Refseq_IDs,c("name"), function(x){
  if(nrow(x)==1){
    return(data.frame(coding_noncoding=as.character(x["coding_noncoding"])))
  }else{
    if(any(x$coding_noncoding == "coding")){
      return(data.frame(coding_noncoding=rep("coding",times=nrow(x))))
    }else{
      return(data.frame(coding_noncoding=rep("noncoding",times=nrow(x))))
    }
  }
})
Refseq_IDs <- temp
Refseq_coding_noncoding <- unique(temp)


######------ write output annotation files -----######
write.table(Refseq_coding_noncoding, file=paste(output,"coding_info.bed",sep=""), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)

