#!/bin/bash
##************************##
## Allelome.PRO2 pipeline ##
##************************##
## Daniel Andergassen
## Tim Hasenbein
## Last modification 04.2025
## Creation: 04.2021


#####*********************** Help message *************************#####
print_help_message(){
    echo "
Usage: bash allelome_pro.sh -i <input_bam> -a <annotation_file> -s <SNP_file> -o <output_directory> [options]

We suggest to give full locations of input files.

Software required:
bedtools (tested on version 2.31.1)
SAMtools (tested on version 1.21)
R (tested on version 4.4.3)
Perl (tested on version 5.22.0)
R package plyr (tested on version 1.8.9)
R package gtools (tested on version 3.9.5)
fetchChromSizes (tested on version 466)
bedToBigBed (tested on version 473)

R packages required:
plyr; gtools
If not installed, the script will try to install them in the default R library path.

Required:
    -i        |    Input sample file (bam format).

    -a        |    Annotation file (BED6 format).

    -s        |    SNP file (BED6 format).

    -o        |    Output directory.

Optional:
    -so       |    Specify if bam file is sorted (1 sorted (default); 0 unsorted).

    -r        |    Min. number of reads for SNPs to be included (default 1).

    -t        |    Total min. number of reads for gene to be included (default 20).

Misc:
    -h        |    Display this help message.
        "
    }

# cmd parsing: flags
while getopts :i:a:s:o:so:r:t:h: flag; do
    case $flag in
        i) sample=$OPTARG;;
        a) annotation=$OPTARG;;
        s) snp_file=$OPTARG;;
        o) outputdir=$OPTARG;;
       so) sorted=$OPTARG;;
        r) minreads=$OPTARG;;
        t) totalreads=$OPTARG;;
    esac
done
# check if required parameters are empty -> if so invoke help
[ -z $sample ] && print_help_message && exit 1
[ -z $annotation ] && print_help_message && exit 1
[ -z $snp_file ] && print_help_message && exit 1
[ -z $outputdir ] && print_help_message && exit 1

# if no value assigned -> default value
sorted="${sorted:-1}"
minreads="${minreads:-1}"
totalreads="${totalreads:-20}"

# check if flag parameters are set
if [ ! -s $sample ]; then
	echo "
Error: Sample bam file ${sample} can not be found."
	exit 0
fi
if [ ! -s $annotation ]; then
	echo "
Error: Annotation file ${annotation} can not be found."
	exit 0
fi
if [ ! -s $snp_file ]; then
	echo "
Error: SNP file ${snp_file} can not be found."
	exit 0
fi
if [ ! -d $outputdir ]; then
	echo "
Error: Output directory ${outputdir} can not be found."
	exit 0
fi

#####********************** Set environment ***********************#####
# get path of the pipe_location
pipe_location=`(dirname "$0")`
# define job directory (delete if already present)
annotation_name=`basename $annotation`
sample_n=`basename "$sample" .bam`
jobdir="${outputdir}/"`date +%Y_%m_%d`"_"$sample_n"_"$annotation_name"_"$minreads
[ -e $jobdir ] && rm -r $jobdir

# define and create directories used throughout script
debug_dir=$jobdir"/debug"
annot=$debug_dir"/annotation"
tmpdir=$jobdir"/sorted_bam_files"
BAM_trim_dir=$debug_dir"/BAM_trim"
BAM_trim_log_sample=$BAM_trim_dir"/trim.log"
BED_dir=$jobdir"/BED_files"
mkdir -p $debug_dir
mkdir -p $BAM_trim_dir
mkdir -p $annot
mkdir -p $BED_dir

#####******************** Start log file *******************#####
# make log file, rm if present
log=$jobdir"/"$sample_n".log"
exec > >(tee -i "$log")
exec 2>&1

# output parameters for control reasons
echo -e "\n\t Start Allelome.PRO run\n"
echo "Input parameters for Allelome.PRO:"
echo "Sample:          $sample"
echo "Annotation:      $annotation"
echo "SNP file:        $snp_file"
echo "Output dir:      $outputdir"
echo "Min reads:       $minreads"
echo "Total min reads: $totalreads"
# check sorting
if [ $sorted == 1 ]; then
  echo 'bamfile sorted: yes
  '
else
  echo 'bamfile sorted: no
  '
fi

# check if BAM file has "chr" prefix
echo -e '---- Checking consistency of 'chr' prefix in the BAM, SNP and ANNOTATION file ----'
if samtools view "$sample" | head -n 1 | cut -f3 | grep -m 1 '^chr'; then
  bam_uses_prefix=TRUE
  echo "Your BAM file has the 'chr' prefix"
  else
  bam_uses_prefix=FALSE
  echo "Your BAM file has no 'chr' prefix"
fi

# function to add 'chr' prefix
add_chr_prefix() {
    local file="$1"
    if ! grep -m 1 '^chr' "$file"; then
        echo "...adding 'chr' prefix to $file"
        # Backup original
        cp "$file" "${file}.bak"
        # Add "chr" to the first column
        awk 'BEGIN{OFS="\t"} {$1 = "chr"$1; print}' "${file}.bak" > "$file"
    else
        echo "$file already has 'chr' prefix"
    fi
}

# function to remove 'chr' prefix
remove_chr_prefix() {
    local file="$1"
    if grep -m 1 '^chr' "$file"; then
        echo "...removing 'chr' prefix from $file"
        # Backup original
        cp "$file" "${file}.bak"
        # Add "chr" to the first column
        awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); print}' "${file}.bak" > "$file"
    else
        echo "$file already does NOT have 'chr' prefix"
    fi
}

# add prefix depending on bam file condition
if [ "$bam_uses_prefix" = TRUE ]; then
    echo "BAM has 'chr' prefix → Ensuring SNP and ANNOTATION match"
    echo "...checking SNP file for 'chr' prefix"
    add_chr_prefix "$snp_file"
    echo "...checking annotation file for 'chr' prefix"
    add_chr_prefix "$annotation"
else
    echo "BAM has no 'chr' prefix → Ensuring SNP and ANNOTATION match"
    echo "...checking SNP file for 'chr' prefix"
    remove_chr_prefix "$snp_file"
    echo "...checking annotation file for 'chr' prefix"
    remove_chr_prefix "$annotation"
fi

#####******************* Executing Allelome.PRO ******************#####
#####=============================================================#####
######------ Interesect SNPs/annotation ------######
# sort bam files if necessary
if [ $sorted = 0 ]; then
	echo "...sorting bam files"
	mkdir -p $tmpdir
  samtools sort -m 6G $sample -o $tmpdir"/"$sample_n"_sorted.bam"
  sample_sorted=$tmpdir"/"$sample_n"_sorted.bam"
  echo -e 'done!\n'
else
	sample_sorted=$sample
fi

snps_in_annot_list=$annot"/snps_overlapping_annotation_list.txt"
snps_in_annot_forjoin=$annot"/snps_overlapping_annotation_forjoin.txt"
annotation_overlapping_snps=$annot"/annotation_overlapping_snps.txt"
annotation_sorted=$annot"/annotation_sorted.bed"
sort -k1,1 -k2,2n $annotation > $annotation_sorted

echo -e '\n---- Intersecting SNP file with annotation file ----'
bedtools intersect  -wa -wb -sorted -a $snp_file -b $annotation_sorted | awk -v OFS="\t" '{split($4,snpdata,""); gsub("chr","",$1); print $1"&"$2,snpdata[1],snpdata[2],$10}' | sort -k 1b,1 -k 4,4 | uniq > $snps_in_annot_forjoin

awk -v OFS="\t" -v list=$snps_in_annot_list -v annot_overlap=$annotation_overlapping_snps '{split($1,pos,"&"); print "chr"pos[1],pos[2] > list; print "chr"pos[1]"\t"$4"\tOK" > annot_overlap; print "chr"pos[1]"\t"pos[2]"\t"pos[2]+1"\t"$2$3"\t0";}' $snps_in_annot_forjoin | sort -k 1,1 -k2,2n | uniq >$annot"/snps_overlapping_annotation.bed"
echo -e '...done!'


######------ Bam Trimming ------######
# start trimming script
echo -e '\n---- Read trimming: ----'
bash $pipe_location"/scripts/bamtrim.sh" $sample_sorted $pipe_location $jobdir $annot"/snps_overlapping_annotation.bed" $snps_in_annot_list
echo -e '...done!'
rm -r $BAM_trim_dir


######------ Pileup file ------######
# Create pileup file: deletion of ^ $ ~ ! and replace chr X/Y chr 20/21, filter by minread threshold set by user, sort file
echo -e '\n---- Creation of pileup file ----'
echo "...creating and formatting pileup file."


## here is the bug fix with gsub(/\^.?/, ""); added -> therefore i remove the ^ASCII flag giving the mapping quality of a SNP if it comes from a read that has the SNP as first base
samtools mpileup -d 50000 -l $snps_in_annot_list $debug_dir"/trimmed_s.bam" | awk -v OFS="\t" -v reads=$minreads '$4 >= reads {gsub("chr","",$1); gsub(/\^.?/, ""); gsub("[\\^|\\$|\\~|\"|!|S|I|D|H|P]","",$5); print $1"&"$2,$4,$5,$6}' | sort -k 1b,1 > $debug_dir"/"$sample_n".pileup"

echo "...joining pileup file with data from SNP list."
join -1 1 -2 1 -o 1.1,2.2,2.3,2.4,1.2,1.3,1.4 $debug_dir"/"$sample_n".pileup" $snps_in_annot_forjoin > $debug_dir"/"$sample_n"_SNP.pileup"
echo -e '...done!'


######------ Allelic scoring ------######
echo -e "\n---- Allelic scoring ----"
Rscript $pipe_location"/scripts/score.R" $jobdir $annotation $debug_dir"/"$sample_n"_SNP.pileup" $totalreads $sample_n $BED_dir
echo -e '...done!'


######------ Creating BigBed file ------######
echo -e "\n---- Creating BigBed file ----"
fetchChromSizes mm10 > $jobdir"/"$sample_n"_mm10.chrom.sizes"
bedToBigBed $BED_dir"/"$sample_n".bed" $jobdir"/"$sample_n"_mm10.chrom.sizes" -type=bed6+3 $BED_dir"/"$sample_n".bb"
rm -f $jobdir"/"$sample_n"_mm10.chrom.sizes"
echo -e '...done!'


######------ Add SNP history to info file ------######
echo -e "\n---- Information about your Allelome.PRO run: ----"
echo -e "Read counts:"
echo -e "Total reads in sample: "`samtools view -c $sample_sorted`
echo -e "Number of reads overlapping at least one SNP in sample: "`samtools view -c $debug_dir"/trimmed_s.bam"`
echo -e "\nSNP history:"
echo -e "Total SNP number: "`cat $snp_file | wc -l`
echo -e "SNP number overlapping with the annotation file: "`cat $snps_in_annot_forjoin | wc -l`
echo -e "SNP number of sample covered by a minimum of $minreads reads or junctions: " `cat $debug_dir"/"$sample_n".pileup" | wc -l`
echo -e "Final SNP number after filtering: "`tail -n +2 $jobdir"/read_count_per_SNP.txt" | cut -f1,2 | sort | uniq | wc -l`
echo -e "\nNumber of loci classified: "`awk 'END { print NR - 1 }' $jobdir"/locus_table.txt"`


echo -e "\n-------------------------------------"
echo -e "     Allelome.PRO run completed"
echo -e "-------------------------------------"
