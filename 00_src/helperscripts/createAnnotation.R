##************************************************##
##  Allelome.PRO / Link generate annotation file  ## 
##************************************************##
## Tim Hasenbein
## Last modification 12.2022
## Creation: 12.2022
## Generate annotation file from RefSeq for Allelome.PRO / Link in BED6 format
## 1. Chromosome, 2. Start Position, 3. End Position, 4. Name, 5 . Coding potential, 6. Strand


######------ set environment -----######
input <- "/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/Y_references/annotation/20180203_RefSeq/"
output <- "/Users/timhasenbein/Desktop/"


######------ get reference data -----######
# gtf
chr <- c(paste(rep("chr",19),1:19,sep=""),"chrX","chrY")
Refseq_gtf <- read.table(file=paste0(input,"20180203_RefSeq.gtf"),sep="\t",header=F,stringsAsFactors=FALSE)
colnames(Refseq_gtf) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
Refseq_gtf <- Refseq_gtf[Refseq_gtf$chr %in% chr,]
# ids
Refseq_IDs <- read.table(file=paste0(input,"20180203_RefSeq_IDs.txt"),sep="\t",header=F,stringsAsFactors=FALSE)
colnames(Refseq_IDs) <- c("transcript_id","chr","strand","start","end","name")
Refseq_IDs <- Refseq_IDs[Refseq_IDs$chr %in% chr,]


######------ replace geneID with gene name -----######
# extract transcript_id from attribute column
Refseq_gtf$transcript_id <- apply(Refseq_gtf,1,function(x){
  return(unlist(strsplit(unlist(strsplit(x["attribute"],split = ";"))[2],split = " transcript_id "))[2])
})
# merge with Refseq_IDs by transcript_id
Refseq_gtf <- merge(Refseq_gtf,Refseq_IDs[,c("transcript_id","name")],by=c("transcript_id"))
Refseq_gtf <- unique(Refseq_gtf)
# replace geneID with name
Refseq_gtf$new_attribute <- apply(Refseq_gtf,1,function(x){
  return(gsub(x["attribute"],pattern = paste("gene_id",x["transcript_id"],sep=" "),replacement = paste("gene_id",x["name"],sep=" ")))
})
Refseq_gtf_output <- Refseq_gtf[,c("chr","source","feature","start","end","score","strand","frame","new_attribute")]


######------ generate annotation file for Allelome.PRO / Allelome.LINK with coding info -----######
Refseq_IDs$score=0
Refseq_IDs=Refseq_IDs[,c("chr","start","end","name","score","strand")]
Refseq_IDs=Refseq_IDs[order(Refseq_IDs$chr,Refseq_IDs$start),]
Refseq_IDs=unique(Refseq_IDs)


######------ generate annotation file for fwd / rev -----######
Refseq_IDs_f <- Refseq_IDs[Refseq_IDs$strand=="+",]
Refseq_IDs_r <- Refseq_IDs[Refseq_IDs$strand=="-",]


######------ write output annotation files -----######
write.table(Refseq_IDs_f, file=paste(output,"annotation_f.bed",sep=""), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)
write.table(Refseq_IDs_r, file=paste(output,"annotation_r.bed",sep=""), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)
write.table(Refseq_IDs, file=paste(output,"annotation_us.bed",sep=""), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)
write.table(Refseq_gtf_output, file=paste(out_folder,"annotation.gtf",sep=""), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)

