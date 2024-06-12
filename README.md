<img width="350" alt="logo" src="https://github.com/AndergassenLab/Allelome.LINK/assets/146940261/587779dd-a554-4c99-b6b9-a2b60e3e1297" align="right">

# Allelome.PRO v2.0     Allelome.LINK

Allelome.PRO is a previously published fully automated bioinformatics pipeline to detect allele-specific expression and histone marks (*Andergassen et. al Nucleic Acids Res. 43, 2015*). Based on heterozygous SNPs, Allelome.PRO assigns sequencing reads to the alleles and classifies your NGS-reads into bi- or monoallelically expressed. By accepting different input data such as RNA-,ChIP-,ATAC- and single-cell-sequencing data, the tool offers a wide range of application. Here, we present Allelome.PRO v2.0, an updated version of the previously published pipeline. Unlike its predecessor, Allelome.PRO v2.0 does not discriminate between ASE loci arising from imprinted or genetic factors. This update streamlines the identification of ASE at the individual level, improving its applicability to diverse biological samples, including human datasets where forward and reverse crosses cannot be obtained. Please find a detailed overview of all updates in the patch notes.

To facilitate the prediction of regulatory interactions and their mode-of-action, we added Allelome.LINK as an extension of Allelome.PRO v2.0. Leveraging the allele-specific information, the pipeline connects ASE loci within user-defined windows in cis and predicts enhancing or repressive effects based on the allelic bias toward identical or opposing alleles. Allelome.LINK offers straightforward execution through a simple one-line command, improving accessibility for diverse users. The output is presented as a tabular list of potential target candidates sorted by linkage score and is accompanied by a BEDPE file for direct visualization, providing an intuitive interface for exploring the results. By default, Allelome.LINK connects all allele-specific loci within the specified window-size. To obtain ncRNA-to-target linkages only, the user has to filter the output file.

#### Publications:
Allelome.PRO v1.0 \
*Andergassen, D., Dotter, C.P., Kulinski, T.M., Guenzl, P.M., Bammer, P.C., Barlow, D.P., Pauler, F.M., and Hudson, Q.J. (2015). Allelome.PRO, a pipeline to define allele-specific genomic features from high-throughput sequencing data. Nucleic Acids Res. 43.*
<br /> \
Allelome.PRO v2.0 / Allelome.LINK \
*Hasenbein T., Andergassen, D.*
<br />

# 1. Installation
To run Allelome.PRO v2.0 / Allelome.LINK no installation is required. The source code is provided as github repository at https://github.com/IPTim/Allelome.LINK. Simply download and unzip the repository and you can start the pipeline via the terminal.

## 1.1 Dependencies 
Software:
```bash
bedtools (≥ version 2.20.1)
SAMtools (≥ version 0.1.19)
R (≥ version 3.1.0)
Perl (≥ version 5.20.0)
fetchChromSizes (≥ version 377)
bedToBigBed (≥ version 377)
```
R packages: 
```bash
plyr; gtools; optparse
# If not installed, Allelome.PRO / Allelome.LINK will try to install them in the default R library path
```

## 1.2 Content of the repository
The github repository contains the following folders and files:

*Table 1* Overview of the repository content
| File  | Description   |
| ------- | ------------- | 
| README.md  | Markdown file with the description of Allelome.PRO v2.0/ Allelome.LINK   | 
| Allelome.PRO | Example run of Allelome.PRO including in-/output files  | 
| Allelome.LINK | Example run of Allelome.LINK including in-/output files  | 
| patch_notes.txt | Release patch notes for Allelome.PRO v2.0  | 
| **src:** | |  
| Allelome.PROv2.0.sh | Main script to run by the user for Allelome.PRO | 
| Allelome.LINK.R   | Main script to run by the user for Allelome.LINK | 
| **src/helperscripts:** | |  
| createAnnotation.sh   | Helperscript to generate annotation file for Allelome.PRO | 
| createCodingfile.sh   | Helperscript to generate coding info file file for Allelome.LINK | 
| createSNPfile_v5.sh   | Helperscript to generate SNP files from mgp.v5.merged.snps_all.dbSNP142.vcf | 
| **src/scripts:** | |
| scripts/bamtrim.pl | Perl script to trim reads covering multiple SNPs to only cover one  | 
| scripts/bamtrim.sh | Bash script to trim reads covering multiple SNPs to only cover one  |   
| scripts/score.R    | R script to generate pileup files and SNP counting. Statisctical scoring of allelic counts and generation of output files | 

<br />

<br />
<img width="250" alt="logo 2" src="https://github.com/AndergassenLab/Allelome.LINK/assets/146940261/6b195065-dd0b-40c2-9eec-2692301e8b46" align="right">


# 3. Run Allelome.PRO v2.0
To run Allelome.PRO v2.0 you need to prepare your input files as described below (3.1). As soon as you have your sample BAM, the annotation and SNP file you can start the pipeline from the command-line by typing:

Usage:
```bash
bash Allelome.PROv2.0.sh -i <input_bam> -a <annotation_file> -s <SNP_file> -o <output_directory> [options]
```

Please make sure that Allelome.PROv2.0.sh is in your PATH or specify the absolute path of the script. For all input files we suggest to give full paths and locations. Besides the required input flags, you have additional options, stated in the table below.

<br />

## 3.1 Input files
Three different files are required as input for Allelome.PRO v2.0:
* Aligned sample file (BAM format)
* Annotation file (BED6 format)
* SNP file (BED6 format)

### 3.1.1 Sample file:
The sample file is a BAM file that stores the aligned NGS-reads from either RNA-, ChIP-, ATAC- or scRNA-sequencing in binary format. The file can be sorted by coordinates or unsorted. If unsorted Allelome.PRO v2.0 can sort the BAM file for you. 

### 3.1.2 Annotation file:
The annotation file is a six column text file containing the posititon information of your loci of interest in BED6 format
(see UCSC format description for more details). The helperscript createAnnotation.sh can be used to generate your own annotation file.
*Table 2* Overview of the BED6 annotation file format
| Column  | Description   | 
| ------- | ------------- | 
| 1   | Chromosome(eg. chr1) | 
| 2   | Start Position  | 
| 3   | End Position | 
| 4   | Name (e.g. gene name)| 
| 5   | Score (not used here)| 
| 6   | Strand (e.g. +,- or . for not defined)| 

### 3.1.3 SNP file:
The SNP file is also a six column text file in BED6 format (see table 2). However, the difference to the annotation file
is that the *name* column consists of two letters representing the two SNPs present at a certain loci. The SNP file is used to determine which SNP is located on which allele. Therefore, the order of the two letter SNPs is important, since the first base indicates the variant at allele 1, while the second base is present at allele 2. The way the SNP file is created defines which allele is "1" and "2". 
SNP positions have to be based on the same reference genome that the BAM file was aligned to (e.g. mm10). One source for mouse SNP data is the FTP site of the Sanger institute (eg. https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/ ). If Allelome.PRO v2.0 is used for hybrid F1 mice samples, the helperscript ./createSNPfile.sh can be used to generate the SNP file via:
```bash
sh create_SNPfile_v5.sh mgp.v5.merged.snps_all.dbSNP142.vcf
```
If Allelome.PRO v2.0 is used with eg. human data, SNPs have to be called and phased in advance to generate the SNP file.
<br />

## 3.2 Input flags for Allelome.PRO.v2.0:
*Table 3.* Overview of the Allelome.PRO v2.0 options.
| Required      | Description   | 
| ------------- | ------------- | 
| -i            | * input sample file (bam format) | 
| -a            | * annotation file (BED6 format)  | 
| -s            | * SNP file (BED6 format)| 
| -o            | * output directory| 
| **Optional:**   |                  | 
| -so           | * specify if bam file is sorted (1 sorted (default); 0 unsorted)| 
| -r            | * min. number of reads to cover a SNP to be included (default 1)| 
| -t            | * min. number of total reads to cover a locus to be included (default 20)| 
| **Misc:**      |                  |
| -h            | * display this help message |

<br />

## 3.3 Output files
Allelome.PRO v2.0 generates a results directory with different files including the allelic bias of your sequencing reads, as well as files for debugging and visualization. The main output files are the **locus_table.txt** file with information on the allelic status of each loci in the annotation and the **<name>.bed** file, to visualize your results in a genome browser. Please, find the detailed information about all output files in Table 3.

*Table 3.*
| Name      | Description   | 
| ------------- | ------------- | 
| locus_table.txt |  Information about the allelic ratio of all informative loci | 
| read_count_per_SNP.txt | Information about the number of reads covering the individual SNPs | 
| <name>.log    | Log file with information about your Allelome.PRO v2.0 run | 
| **./BED_files:**  |                  | 
| <name>.bed    | BED6-file containing all loci, color-coded according to their allelic ratio| 
| <name>.bb     | Big bed file of <name>.bed | 
| **/debug**      | Folder containing files created during the run for debugging |
| <name>_SNP.pileup    |  | 
| <name>.pileup     |  | 
| trimmed_s.bam     |  | 
| **/debug/annotation**      | Folder containing files created during the run for debugging |
| annotation_overlapping_snps.txt    |  | 
| annotation_sorted.bed  |  | 
| snps_overlapping_annotation_forjoin.txt  |  | 
| snps_overlapping_annotation_list.txt  |  | 
| snps_overlapping_annotation.bed  |  | 

<br />
<br />
<img width="250" alt="logo 2" src="https://github.com/AndergassenLab/Allelome.LINK/assets/146940261/e02b66ad-6bc4-410e-81ef-5d83b355b98b" align="right">

# 4. Run Allelome.LINK
To run Allelome.LINK you just need your *locus_table.txt* file from the Allelome.PRO v2.0 run. To start Allelome.LINK, type in your command line:

Usage: 
```bash
Rscript Allelome.LINK.R -i <input_locus_table.txt> -o <output_directory> [options] 
```
Again, please make sure that Allelome.LINK is in your PATH or give the full direction to the file. Allelome.LINK is linking loci within a given window size in *cis*. You can specifiy your the genomic window by the --window-size or -w flag which is the number of basepairs up- and downstream of the locus (see table 4). The default value is 100Kb in each direction from a given locus. Further you can specify the number of total reads that have to cover a gene to consider it for the analysis (--total-reads/-r; default: 20), as well as the --allelelic-bias/-b to define the cut-off that distinguishes bialleic from allele-specific expression. The default value is 0.7, which means that at least 70% of the reads from a genes have to be expressed from one allele to consider it as imbalanced. Please find all options described in Table 4.

<br />

## 4.1 Input files
* Locus table (.txt format) as given by Allelome.PRO.
To run Allelome.LINK you simply need the **locus_table.txt** file from your previous Allelome.PRO run.

<br />

## 4.2 Input flags for Allelome.LINK
*Table 4*
| Required      | Description   | 
| ------------- | ------------- | 
| -input         -i   | * input locus_table.txt (as given by Allelome.PRO v2.0) | 
| **Optional:**   |                  | 
| --name         -n   | * sample name (default date and time)| 
| --window-size  -w | * window range to draw links (in Kb; default ±100)| 
| --total-reads  -r | * total read cut-off to consider allele-specific genes (default 20)| 
| --allelic-bias -b| * cut-off to define genes with an allelic bias (default 0.7)| 
| --duplicates -d  | * remove mirrored duplicates (default TRUE) |
| --output       -o | * output directory (default ./) | 
| **Misc:**      |                  |
| --help  -h   | * display this help message |

<br />

## 4.3 Output files 
Allelome.LINK generates a results directory including different output files. The main output files are the **links_table.txt** with information on the linked loci and the **.bedpe** file, to visualize the links in a genome browser. Please, find the detailed information about all output files in Table 5.

*Table 5.*
| Name      | Description   | 
| ------------- | ------------- | 
| <name>_links_table.txt | Information about the predicted regulatory association between nearby loci. | 
| <name>_links_full_table.txt  | As <name>_links_table.txt but including the read and allele-specific information |                                  | 
| <name>.log    | Log file with information about your Allelome.LINK run | 
| **/BED_files**  |                  | 
| <name>.bed | BED6-file containing all loci passed, color-coded according to their allelic bias | 
| <name>.bedpe | BEDPE file for link representation, color-coded accroding to mechanism | 
| <name>_repressing.bedpe | BEDPE file for IGV browser visualization, repressive links only | 
| <name>_enhancing.bedpe | BEDPE file for IGV browser visualization, enhancing links only | 
<br />

# 5. Resource generated by <cite>
To access the resource generated by <cite> download json files for mice and human samples from the /resource/ folder. 
Candidate predictions can then be viewed via the Integrative Genomics Viewer (IGV) at https://igv.org/app/. 
Please load json files as session. Upon loading the json file of interest you can see the following output:

## 5.1 Example mice: Heart 
<img width="951" alt="mice_example" src="https://github.com/AndergassenLab/Allelome.LINK/assets/146940261/d79d4879-a786-48ea-8294-1dc6f06e08fb">


IGV: James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011). A public access version is also available: PMC3346182.


## Contact
For any question, issue, bug or recommendation, please submit a github issue.
