##*******************##
## Allelome.LINK run ##
##*******************##
## Project: Allelome.Link
## Tim Hasenbein
## Last modification 06.2024
## Creation: 06.2024
## Example run for Allelome.LINK


######------ Set environment ------######
pipeline="./00_src/Allelome.LINK.R"
input="./02_Allelome.PRO_v2.0/2024_06_05_sample_br_9w_BxC_annotation.bed_1/locus_table.txt"
total_read=20
bias=0.7
window=100
output_dir="./03_Allelome.LINK"

######------ Run Allelome.PRO v2.0 ------######
Rscript $pipeline -i $input -r $total_read -b $bias -w $window -o $output_dir