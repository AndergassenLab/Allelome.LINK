#!/usr/bin/env Rscript
##*************************##
## Allelome.LINK pipeline  ## 
##*************************##
## Tim Hasenbein
## Last modification 07.2024
## Creation: 02.2022
## Allelome.Link pipeline


#####******************** Command line parsing *******************#####
argv <- commandArgs(trailingOnly = FALSE)


#####*********************** Help message *************************#####
help_message <- 
  "\n
*********************** Allelome.Link *****************************
Usage: Rscript Allelome.LINK.R -i <input_locus_table.txt> -o <output_directory> [options]\n 
R packages required:
optparse
If not installed, Allelome.LINK will try to install them in the default R library path.\n
Required:
    --input         | -i    Input locus table (as given by Allelome.PRO).
    
Optional:
    --name          | -n    Sample name (default date and time).
    
    --window-size   | -w    Window range to draw links (in Kb; default 500).
    
    --total-reads   | -r    Total read cut-off to consider biased genes (default 50).
    
    --allelic-bias  | -b    Cut-off to define genes with an allelic bias (default 0.7).
    
    --duplicates    | -d    Remove mirrored duplicates (default TRUE).
    
    --output        | -o    Output directory (default ./).
  
Misc:
    --help          | -h    Display this help message.\n\n"
myhelp <- (sum(c("--help", "-h") %in% argv) >= 1) 
if (myhelp) {
  cat(help_message)
  quit()
}
# check if requirements are available
req <- (sum(c("--input", "-i") %in% argv) < 1) # add name
if (req) {
  cat(help_message)
  quit()
}


#####************************* Packages ***************************#####
packages <- function(requirements,quiet=FALSE){
  has <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=c(1:2))
    install.packages(requirements[!has],repos = "http://cran.us.r-project.org")
  }
  if(quiet){
    for(r in requirements){suppressMessages(require(r,character.only=TRUE))}
  }else for(r in requirements){message(paste(r,suppressMessages(require(r,character.only=TRUE)),sep=': '))}
}
packages(c('optparse')) 


#####********************* Pipeline options ***********************#####
option_list <- list(
  make_option(c("--input","-i"),type = "character",dest = "input",default=''),
  make_option(c("--name","-n"),type = "character",dest = "name",default=paste0("AL_run_",Sys.time())), 
  make_option(c("--window-size","-w"),type = "numeric",dest = "window",default=500),
  make_option(c("--total-reads","-r"),type = "numeric",dest = "read",default=50),
  make_option(c("--allelic-bias","-b"),type = "numeric",dest = "bias",default=0.7),
  make_option(c("--duplicates","-d"),type = "logical",dest = "dups",default=T),
  make_option(c("--output","-o"),type = "character",dest = "out",default='./')
)
# If help is not invoked, start processing the input
opt <- parse_args(OptionParser(option_list=option_list))


#####********************** Set environment ***********************#####
output <- paste0(opt$out,"/",opt$name,"/")
window <- opt$window*1000 
bias_up <- opt$bias
bias_down <- 1 - opt$bias
dir.create(output) 


#####******************** Start log file *******************#####
log <- file(paste0(output,"log.txt"), open = "wt")
sink(file = log,type="output")
sink(file = log,type="message")


#####******************* Executing Allelome.Link ******************#####
#####==============================================================#####
message("\n------------  Starting Allelome.Link:  ------------ ")
message(paste0("Input: ",opt$input))
message(paste0("Name: ",opt$name))
message(paste0("Allelic ratio: ",bias_up,"/",bias_down))
message(paste0("Window-size: ",window/1000,"Kb")) 
message(paste0("Total min. read: ",opt$read))
message(paste0("Output: ",opt$out))
message("--------------------------------------------------- ")


######------ Get input ------###### 
ap <- read.table(paste0(opt$input), sep='\t', row.names = NULL, header=T)
number_of_input_loci <- nrow(ap)


######------ Get biased loci -----###### 
message("\n---- Get biased loci ----")
ap_wip <- ap[ap['total_reads'] >= opt$read & (ap['chr'] != 'chrX') & (ap['chr'] != 'chrY'),] 
number_of_loci_minread <- nrow(ap_wip)
# filter biased elements
ap_wip <- ap_wip[ap_wip['allelic_ratio'] >= bias_up | ap_wip['allelic_ratio'] <= bias_down,]
number_of_biased_loci <- nrow(ap_wip)
message("done.")


######------ Link loci within windows -----######  
message("\n---- Generate links ----")
linkages <- data.frame() 
for(row in 1:nrow(ap_wip)){
  # get base gene dataframe
  base <- ap_wip[row,]
  # get window size
  downstream <- base[,'start'] - window
  upstream <- base[,'end'] + window
  # filter couples and create results
  matrix <- ap_wip[base[,'chr'] == ap_wip[,'chr'],] 
  matrix <- matrix[base[,'name'] != matrix[,'name'],] 
  couple <- matrix[(matrix[,'start'] >= downstream & matrix[,'start'] <= upstream) | 
                     (matrix[,'end'] >= downstream & matrix[,'end'] <= upstream) |
                     (matrix[,'end'] >= downstream & matrix[,'start'] <= upstream),] 
  colnames(base) <- c("chr_base", "start_base", "end_base", 
                      "name_base", "A1_reads_base","A2_reads_base",
                      "total_reads_base","allelic_score_base","allelic_ratio_base")   
  colnames(couple) <- c("chr_target", "start_target", "end_target", 
                        "name_target", "A1_reads_target","A2_reads_target",
                        "total_reads_target","allelic_score_target","allelic_ratio_target")   
  ifelse(nrow(couple) == 0, link <- NA, link <- cbind(base, couple, row.names = NULL)) 
  ifelse(is.na(link), linkages <- linkages, linkages <- rbind(linkages, link, row.names = NULL))
} 
message("done.")


######------ Assign mechanism -----###### 
if(nrow(linkages)!=0) { 
  message("\n---- Assign mechanism ----")
  linkages$mechanism <- NA 
  linkages$mechanism[(linkages$allelic_ratio_base>=bias_up & linkages$allelic_ratio_target>=bias_up)|
                     (linkages$allelic_ratio_base<=bias_down & linkages$allelic_ratio_target<=bias_down)] <- "enhancing"
  linkages$mechanism[is.na(linkages$mechanism)] <- "repressing"
  message("done.")


  ######------ Add linkage score -----######
  score <- data.frame("base"=abs(linkages$allelic_score_base),"target"=abs(linkages$allelic_score_target))
  score$min <- apply(score, 1, min)
  ratio <- data.frame("base"=abs(linkages$allelic_ratio_base - 0.5), "target"=abs(linkages$allelic_ratio_target - 0.5)) 
  ratio$delta <- 1-(abs(ratio$base - ratio$target))
  linkages$linkage_score <- round(x=log10(score$min+1) * ratio$delta,digits = 2) 


  ######------ Remove mirrored duplicates -----###### 
  if(opt$dups){
    linkages$grp <- paste(pmax(linkages$name_base, linkages$name_target), pmin(linkages$name_base, linkages$name_target), sep = "_")   
    linkages <- linkages[!duplicated(linkages$grp),]
    linkages <- linkages[,c(1:20)]
  }


  ######------ Create bedpe file ------###### 
  pos1 <- round((linkages$start_base + linkages$end_base)/2, digits=0)
  pos2 <- round((linkages$start_target + linkages$end_target)/2, digits=0)
  bedpe <- data.frame(V1 = linkages$chr_base,V2 = pos1,V3 = pos1,V4 = linkages$chr_target,V5 = pos2,V6 = pos2,V7 = paste0(linkages$name_base,"_",linkages$name_target),V8 = linkages$linkage_score,V9 = ".",V10 = ".",V11 = linkages$mechanism)
  # add color according to mechanism
  bedpe[bedpe$V11=="enhancing","V11"] <- "71,133,119"
  bedpe[bedpe$V11=="repressing","V11"] <- "188,34,13"
  # workaround: add this line to bedpe to accept colors
  col_header <- data.frame(V1="chr1",V2="x1",V3="x2",V4="chrom2",V5="start2",V6="end2",V7="name",V8="score",V9="strand1",V10="strand2",V11="color")
  bedpe <- rbind(col_header,bedpe)
  colnames(bedpe) <- c("track='interact'","thickness=3","","","","","","","","","")
  write.table(bedpe,paste0(output,opt$name,".bedpe"),sep="\t",quote=F,row.names = F, col.names = T)
  # filter bedpe for just enhancing/repressing
  bedpe_enhancing <- bedpe[bedpe[,11]=="71,133,119"|bedpe[,11]=="color",] 
  bedpe_repressing <- bedpe[bedpe[,11]=="188,34,13"|bedpe[,11]=="color",] 
  write.table(bedpe_enhancing,paste0(output,opt$name,"_enhancing.bedpe"),sep="\t",quote=F,row.names = F, col.names = T)
  write.table(bedpe_repressing,paste0(output,opt$name,"_repressing.bedpe"),sep="\t",quote=F,row.names = F, col.names = T)


  ######------ Create bed file ------###### 
  bed <- ap[c(1:3,7,9)]
  bed$name <- paste0(ap$name,"_reads:",ap$A1_reads,"/",ap$A2_reads,"_score:",ap$allelic_score,"_ratio:",ap$allelic_ratio)
  bed$score <- "1000"
  bed$strand <- "."
  bed$thickStart <-  paste0(bed$start)
  bed$thickEnd <- paste0(bed$end)
  bed$itemRgb[bed$allelic_ratio >= 0.6 & 
              bed$allelic_ratio < 0.8] <- "207,136,131"
  bed$itemRgb[bed$allelic_ratio >= 0.8 & 
              bed$allelic_ratio <= 1] <- "139,25,19"
  bed$itemRgb[bed$allelic_ratio <= 0.4 & 
              bed$allelic_ratio > 0.2] <- "140,164,202"
  bed$itemRgb[bed$allelic_ratio <= 0.2 & 
              bed$allelic_ratio >=0] <- "48,74,153"
  bed$itemRgb[bed$allelic_ratio > 0.4 & 
              bed$allelic_ratio < 0.6] <- "28,100,45"
  bed$itemRgb[bed$total_reads < opt$read] <- "80,80,80"
  bed <- bed[,c(1:3,6:11)]
  colnames(bed) <- c("track",paste0("name=",opt$name),paste0("description=",opt$name),"visibility=2","itemRgb=On","useScore=1","","","")
  write.table(bed, paste0(output,opt$name,".bed"), quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)


  ######------ Write output -----###### 
  message("\n---- Write output ----")
  linkages <- linkages[order(-linkages$linkage_score),]
  write.table(linkages, paste0(output,opt$name,"_links_full_table.txt"),sep="\t",quote=F,row.names = F) 
  write.table(linkages[,c("chr_base","start_base","end_base","name_base","name_target","start_target","end_target","linkage_score")], paste0(output,opt$name,"_links_table.txt"),sep="\t",quote=F,row.names = F) 
  message("done.")


  ######------ Summary -----###### 
  message("\n---- Information about your Allelome.Link run: ----")
  message("Number of informative loci:            ",number_of_input_loci) 
  message("Number of loci passing read cut-off:   ",number_of_loci_minread) 
  message("Number of loci with allelic imbalance: ",number_of_biased_loci)
  message("Number of unique links:                ",nrow(linkages)) 
  message("   Thereof enhancing:                  ",sum(linkages$mechanism=="enhancing"))
  message("   Thereof repressing:                 ",sum(linkages$mechanism=="repressing"),"\n")
  message("\n---- Top 10 links within your data: ----")
  print(linkages[1:10,c("name_base","name_target","linkage_score")]) 
  message("\n--------------------------------------")
  message("      Allelome.Link run completed") 
  message("--------------------------------------\n")
  sink()
} else {
  message("\nAllelome.Link run completed: No linkages found.")
}






