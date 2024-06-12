##******************##
## Allelome.PRO run ##
##******************##
## Project: Allelome.Link
## Tim Hasenbein
## Last modification 06.2024
## Creation: 06.2024
## Example run for Allelome.PRO v2.0


######------ Set environment ------######
pipeline="./00_src/Allelome.PROv2.0.sh"
input="./01_input/sample_br_9w_BxC.bam"
annotation="./01_input/annotation.bed"
snps="./01_input/SNPfile_BxC.bed"
min_read=1
total_read=1
output_dir="./02_Allelome.PRO_v2.0/"


######------ Run Allelome.PRO v2.0 ------######
bash $pipeline -i $input -a $annotation -s $snps -r $min_read -t $total_read -o $output_dir


