##************************##
## Allelome.PRO2/LINK run ##
##************************##
## Project: Allelome.LINK
## Tim Hasenbein
## Last modification 05.2025
## Creation: 05.2025
## Performe an analysis run on a subset data set using Allelome.PRO2 and Allelome.LINK


######------ Create environment ------######
conda env create -f ./01_test_folder/00_environment/environment.yml --name AllelomePRO2_test
conda activate AllelomePRO2_test


######------ Run Allelome.PRO2 ------######
bash "./00_src/Allelome.PRO2.sh" \
  -i ./01_test_folder/01_input/01_BAM_chr17.bam \
  -s ./01_test_folder/01_input/02_SNPS_chr17.bed \
  -a ./01_test_folder/01_input/03_ANNOTATION_chr17.bed \
  -r 1 \
  -t 20 \
  -o ./01_test_folder/02_output


######------ Capture name of Allelome.PRO2 output ------######
# Capture the newest output folder name (assuming it's the most recently created/modified)
latest_output=$(ls -td ./01_test_folder/02_output/*/locus_table.txt | head -n 1)


######------ Run Allelome.LINK ------######
Rscript "./00_src/Allelome.LINK.R" \
  -i $latest_output \
  -r 20 \
  -b 0.7 \
  -w 100 \
  -n test_run \
  -o ./01_test_folder/02_output
