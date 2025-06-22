##*********************************##
## Scoring script for Allelome.PRO ##
##*********************************##
## Project: Allelome.PRO
## Daniel Andergassen
## Tim Hasenbein
## Last modification 02.2024
## Creation: 04.2021


######------ Set environment ------###### 
packages <- function(requirements,quiet=FALSE){
  has <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=c(1:5))
    install.packages(requirements[!has], repos="https://cran.uni-muenster.de/", dependencies = T)
  }
  if(quiet){
    for(r in requirements){suppressMessages(require(r,character.only=TRUE))}
  }else for(r in requirements){message(paste(r,suppressMessages(require(r,character.only=TRUE)),sep=': '))}
}
packages(c('plyr','gtools'))
path <- commandArgs(TRUE) 
job_path <- path[1] 
annotation_path <- path[2] 
pileup <- path[3] 
totalreads <- as.numeric(path[4])
name <- path[5]
bed_dir <- path[6]


######------ Filter pileup and read count ------###### 
pileup <- read.delim(pileup, header = F, stringsAsFactors = FALSE, sep = " ")
# remove indels
num_insertions <- nrow(pileup[grepl('\\+', pileup$V6), ]) 
num_deletions <- nrow(pileup[grepl('\\-', pileup$V6), ]) 
num_total <- nrow(pileup) 
pileup <- pileup[!grepl('\\d+', pileup$V6), ] 
message(paste0("\n...removed ",num_insertions," (",round((num_insertions/num_total)*100, 3),"%) insertions and ",num_deletions," (",round((num_deletions/num_total)*100,3),"%) deletions.")) 
# double check if number of SNVs and QS overlap
pileup$snvs <- nchar(pileup$V6) 
pileup$asc <- nchar(pileup$V7)
pileup <- pileup[pileup$snvs == pileup$asc, ] 
pileup <- cbind(pileup, data.frame(t(sapply(strsplit(pileup$V1,'&'), `[`)))) 
pileup <- pileup[,c("X1","X2","V2","V3","V4","V5","V6","V7")]
colnames(pileup) <- c("chr", "pos", "A1", "A2", "name", "reads_total", "SNPs", "phred")
pileup$chr <- paste0('chr',pileup$chr)
# count A1, A2 reads and update total reads
message("...filter for quality and counting reads over each locus.")
read.count.function <- function(x) {
  base <- strsplit(as.character(x$SNPs), split="")
  qual <- asc(as.character(x$phred), simplify = F)
  fil <- do.call(rbind, Map(data.frame, allele=base,quality=qual))
  fil <- fil[fil$allele != ">" & fil$allele != "<" & fil$allele != "*", ] 
  x$SNPs <- toupper(paste(fil$allele, collapse =""))
  return(cbind(lengths(regmatches(x$SNPs, gregexpr(x$A1, x$SNPs))), lengths(regmatches(x$SNPs, gregexpr(x$A2, x$SNPs)))))
}
read_count <- ddply(pileup,c("chr", "pos", "A1", "A2", "name"), read.count.function)
colnames(read_count) <- c("chr", "pos", "A1", "A2", "name", "A1_reads", "A2_reads")
read_count <- read_count[read_count$A1_reads + read_count$A2_reads != "0", ] 
write.table(read_count, paste(job_path,"/read_count_per_SNP.txt",sep=""), quote=FALSE, sep = "\t", col.names = TRUE, row.names =FALSE)
rm(pileup)


######------ Scoring script ------###### 
# allelic score funtion
calc.allelic.score <- function(x) {
  if (all(x==0) || any(is.na(x)))
    return(0)
  if (sum(x) < 10)
    return(0)
  allelic.score <- log10(pbinom(min(matrix(x,ncol=2)), x[1]+x[2], 0.5, lower.tail=TRUE, log = FALSE)) 
  ifelse(x[1] > x[2], -allelic.score,allelic.score) 
}


######------  Datamatrix ------###### 
message("...processing datamatrix.")
read_count$pos <- as.numeric(as.character(read_count$pos)) 
annotation <- read.delim(annotation_path,header=FALSE,stringsAsFactors=FALSE)
annotation <- unique(annotation)
colnames(annotation) <- c("chr","start","end","name","bedscore","strand") 
# makes longest possible transcript (locus) out of all the isoforms of a gene 
annotation=ddply(annotation, c("chr","strand","name"), function(x){ 
  return(cbind(min(x[,"start"])[1],max(x[,"end"])[1]))
})
colnames(annotation) <- c("chr","strand","name","start","end")
annotation$start <- as.integer(annotation$start)
annotation$end <- as.integer(annotation$end)
# create data_matrix data frame 
data_matrix <- merge(read_count,annotation,by=c("chr","name"))
data_matrix <- unique(data_matrix)
data_matrix <- data_matrix[,c(1:2,8:10,3:7)] 
data_matrix <- data_matrix[(data_matrix$pos <= data_matrix$end & data_matrix$pos >= data_matrix$start),]
data_matrix <- data_matrix[order(data_matrix$chr,data_matrix$pos),]
rm(read_count)


######------  Locus table ------###### 
# summing up all reads over a SNP within a gene
message("...summing up all reads over a SNP within a gene.")
locus_table <- ddply(data_matrix, c("chr","start","end","name"),function(x) data.frame(A1_reads=sum(x[9],na.rm=T),A2_reads=sum(x[10],na.rm=T)))
# ordering data by chr and start of genes
locus_table <- locus_table[order(locus_table$chr,locus_table$start), ]
locus_table <- unique(locus_table)
# total counts per gene 
locus_table$total_reads <- locus_table$A1_reads + locus_table$A2_reads


######------  Calculating allelic score ------###### 
# Loci: Calculate Imprinting Score
message("...calculating allelic score.")
locus_table$allelic_score <- round(apply(locus_table[,c("A1_reads","A2_reads")],1,calc.allelic.score),3)
# Loci: Replacing +/-infinitive with the highest/lowest finite value
max_value <- range(abs(locus_table$allelic_score),finite=TRUE)[2]
locus_table$allelic_score[locus_table$allelic_score == "Inf"] <- max_value
locus_table$allelic_score[locus_table$allelic_score == "-Inf"] <- -max_value


######------  Calculating allelic ratio ------###### 
# Loci: Calculates the allelic ratio
message("...calculating the allelic ratio.")
locus_table$allelic_ratio <- round(apply(locus_table[,c("A1_reads","total_reads")],1,function(x) {return(x["A1_reads"]/x["total_reads"])}),2)


######------  Apply total reads cut-off and write locus table ------###### 
locus_table_filtered <- locus_table[locus_table$total_reads >= totalreads,]
message("...writing output.")
write.table(locus_table_filtered, paste(job_path,"/locus_table.txt",sep=""), quote=FALSE, sep = "\t", col.names = TRUE, row.names =FALSE)
message("...done!")


######------  Creating bed file ------###### 
message("\n---- Creating bed file ----")
bedfile <- locus_table[c(1:3,7,9)]
bedfile$name <- paste0(locus_table$name,"_reads:",locus_table$A1_reads,"/",locus_table$A2_reads,"_score:",locus_table$allelic_score,"_ratio:",locus_table$allelic_ratio)
bedfile$score <- "1000"
bedfile$strand <- "."
bedfile$thickStart <-  paste0(locus_table$start)
bedfile$thickEnd <- paste0(locus_table$end)
# assign color
bedfile$itemRgb[bedfile$allelic_ratio >= 0.7 & 
                  bedfile$allelic_ratio < 0.8] <- "207,136,131"
bedfile$itemRgb[bedfile$allelic_ratio >= 0.7 & 
                  bedfile$allelic_ratio <= 1] <- "139,25,19"
bedfile$itemRgb[bedfile$allelic_ratio <= 0.3 & 
                  bedfile$allelic_ratio > 0.2] <- "140,164,202"
bedfile$itemRgb[bedfile$allelic_ratio <= 0.2 & 
                  bedfile$allelic_ratio >=0] <- "48,74,153"
bedfile$itemRgb[bedfile$allelic_ratio > 0.3 & 
                  bedfile$allelic_ratio < 0.7] <- "28,100,45"
bedfile$itemRgb[bedfile$total_reads < totalreads] <- "80,80,80"
bedfile <- bedfile[,c(1:3,6:11)]
colnames(bedfile) <- c("track",paste0("name=",name),paste0("description=",name),"visibility=2","itemRgb=On","useScore=1","","","")
write.table(bedfile, paste0(bed_dir,"/",name,".bed"), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)

