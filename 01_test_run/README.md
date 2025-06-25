## Test end-to-end run for Allelome.PRO2 and Allelome.LINK

Simply download the repository and get started:

1. Create environment: Create an environment with all necessary software packages.
2. Run Allelome.PRO2: Executes the Allelome.PRO2 pipeline on the test data set.
3. Run Allelome.LINK Executes the Allelome.LINK pipeline on the test data set.

To execute give the full path to the files or execute from main.

## Sample description
- BAM file: RNA-seq data (unstranded) of the heart from 9-week old F1 hybrid mice (BL6xCAST cross).
- SNP file: Generated for the BL6xCAST cross.
- Annotation: RefSeq annotation.
All files are subset for chromosome 17.

## 1. Create environment
```bash
######------ 1. Create environment ------######
# Does not work for arm64
conda env create -f ./01_test_run/00_environment/environment.yml --name AllelomePRO2_test
conda activate AllelomePRO2_test
```

## 2. Run Allelome.PRO2
```bash
######------ 2. Run Allelome.PRO2 ------######
bash "./00_src/Allelome.PRO2.sh" \
  -i ./01_test_run/01_input/01_BAM_chr17.bam \
  -s ./01_test_run/01_input/02_SNPS_chr17.bed \
  -a ./01_test_run/01_input/03_ANNOTATION_chr17.bed \
  -r 1 \
  -t 20 \
  -o ./01_test_run/
```

## 3. Run Allelome.LINK
```bash
######------ Capture name of Allelome.PRO2 output ------######
# Capture the newest output folder name (assuming it's the most recently created/modified)
latest_output=$(ls -td ./01_test_run/*/locus_table.txt | head -n 1)

######------ 3. Run Allelome.LINK ------######
Rscript "./00_src/Allelome.LINK.R" \
  -i $latest_output \
  -r 20 \
  -b 0.7 \
  -w 100 \
  -n AllelomeLINK \
  -o ./01_test_run/
```

### Subsequently, you can load the .bed and .bedpe files into the IGV genome browser

<img width="1029" alt="Screenshot 2025-06-25 at 14 15 24" src="https://github.com/user-attachments/assets/0dd63112-7c2e-49e2-9a52-d6dd1d8e508f" />




