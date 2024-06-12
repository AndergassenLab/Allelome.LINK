#!/usr/bin/env Rscript
##*****************************************************************##
## Allelome.LINK pipeline: Filter ncRNA to pcRNA interactions only ## 
##*****************************************************************##
## Tim Hasenbein, Daniel Andergassen
## Last modification 12.2022
## Creation: 12.2022
## This script takes the Allelome.LINK output as input and 
## generates new Allelome.LINK output files including the nc-pc links only
## Input requirements: 1. Alllelome.LINK output, 2. Coding info file
## For RNAseq data only


#####******************** Command line parsing *******************#####
argv <- commandArgs(trailingOnly = FALSE)


#####*********************** Help message *************************#####
help_message <- 
  "\n
*********************** Allelome.Link helperscript *****************************
Helperscript for Allelome.LINK to filter the output for noncoding to protein-coding links only.\n 
Usage: Rscript filter_ncpc_links.R -i <Allelome.LINK output> -c <coding info> -o <output_directory> \n 
R packages required:
optparse
If not installed, Allelome.LINK will try to install them in the default R library path.\n
Required:
    --input         | -i    Input links full table (as given by Allelome.LINK).
    
    --bedpe         | -b    Input bedpe file (as given by Allelome.LINK).
    
    --coding        | -c    Info file about coding potential of loci.
    
Optional:
    --name          | -n    Name of output files (default ncpc).
    
    --output        | -o    Output directory (default ./).
  
Misc:
    --help          | -h    Display this help message.\n\n"
myhelp <- (sum(c("--help", "-h") %in% argv) >= 1) 
if (myhelp) {
  cat(help_message)
  quit()
}
req <- (sum(c("--input", "-i","--bedpe", "-b","--coding", "-c") %in% argv) < 3)
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
  make_option(c("--bedpe","-b"),type = "character",dest = "bedpe",default=''),
  make_option(c("--coding","-c"),type = "character",dest = "coding",default=''), 
  make_option(c("--name","-n"),type = "character",dest = "name",default='ncpc'), 
  make_option(c("--output","-o"),type = "character",dest = "out",default='./')
)
opt <- parse_args(OptionParser(option_list=option_list))


#####******************* Executing script ******************#####
#####=======================================================#####
# --- dev --- #
#input <- "/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/03_project_allelome/03_bodymap_9w/03_allelomelink_9w/AL_r20_b0.6_w500/He_merged_RNA_m20_b0.6_w500_test_2/He_merged_RNA_m20_b0.6_w500_test_2_links_full_table.txt"
#bedpe <- "/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/03_project_allelome/03_bodymap_9w/03_allelomelink_9w/AL_r20_b0.6_w500/He_merged_RNA_m20_b0.6_w500_test_2/He_merged_RNA_m20_b0.6_w500_test_2.bedpe"
#coding_info <- "/Users/timhasenbein/Desktop/coding_info.bed"
#al_tab <- read.table(input, header=T)
#bedpe <- read.table("/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/03_project_allelome/03_bodymap_9w/03_allelomelink_9w/AL_r20_b0.6_w500/He_merged_RNA_m20_b0.6_w500_test_2/He_merged_RNA_m20_b0.6_w500_test_2.bedpe", header=F, fill = TRUE)
#bedpe <- bedpe[-c(1,2),]
#coding_info <- read.table(coding_info, header=F)
#-------------#
######------ get input files -----###### 
al_tab <- read.table(opt$input, header=T)
bedpe <- read.table(opt$bedpe, header=F, fill = TRUE)
bedpe <- bedpe[-c(1,2),]
coding_info <- read.table(opt$coding, header=F)


######------ add coding info -----###### 
colnames(coding_info) <- c("name_base","coding_info_base")
al_tab <- merge(al_tab,coding_info,by="name_base",all.x = T)
colnames(coding_info) <- c("name_target","coding_info_target")
al_tab <- merge(al_tab,coding_info,by="name_target",all.x = T)


######------ add type of linkage-----###### 
al_tab$interaction <- NA
al_tab$interaction[al_tab$coding_info_base=="noncoding" & 
                   al_tab$coding_info_target=="noncoding"] <- "nc_nc"
al_tab$interaction[al_tab$coding_info_base=="coding" & 
                   al_tab$coding_info_target=="coding"] <- "pc_pc"
al_tab$interaction[is.na(al_tab$interaction)] <- "nc_pc"


######------ filter for nc_pc links from linkage table -----###### 
nc_pc <- al_tab[al_tab$interaction == "nc_pc",]
nc_pc$link <- paste0(nc_pc$name_base,"_",nc_pc$name_target)


######------ filter for nc_pc links from bedpe file -----###### 
bedpe_ncpc <- bedpe[bedpe$V7 %in% nc_pc$link,]
colnames(bedpe_ncpc) <- c("track='interact'","thickness=3","","","","","","","","","")
bedpe_ncpc <- rbind(c("chr1","x1","x2","chrom2","start2","end2","name","score","strand1","strand2","color"), bedpe_ncpc)


######------ write output -----###### 
write.table(nc_pc[,c(1:20)],paste0(opt$out,"/",opt$name,"_links_full_table.txt"),quote=F,col.names = T, row.names = F, sep = "\t")
write.table(bedpe_ncpc,paste0(opt$out,"/",opt$name,".bedpe"),quote=F,col.names = T, row.names = F, sep = "\t")



