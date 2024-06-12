#!/bin/bash
set -e
# additional software requirements:
# bedtools>=2.20.1
# samtools>=0.1.19

bam=$1
pipeline_location=$2
job_dir=$3
snps_overlapping_annotation=$4
snps_in_annot_list=$5
scripts_dir=$pipeline_location"/scripts"
debug_dir=$job_dir"/debug"
annot=$debug_dir"/annotation"
BAM_trim_dir=$debug_dir"/BAM_trim"

echo -e "Parameter information:"
echo "sample: $bam"
echo "pipeline location: $pipeline_location"
echo "job directory: $job_dir"
echo "SNPs overlapping annotation file: $snps_overlapping_annotation"
echo "SNPs overlapping annotation list: $snps_in_annot_list"
echo "debug directory: $debug_dir"
echo "annotation directory: $annot"
echo "tmp bamtrim directory: $BAM_trim_dir"
echo "script directory: $scripts_dir"

echo -e "\nStart trimming process..."
echo -e "...intersect SNPs overlapping annotation with reads overlapping SNPs. "
# ensure that the snp bedfile is sorted the same way as the bam file; bam sort order is determined by the bam header
samtools view -H $bam | grep SN: | awk -v OFS="\t" '{split($2,n,":"); split($3,s,":"); print n[2],s[2]}' > $BAM_trim_dir/"chr_order_in_bam_header.txt"
cat $snps_overlapping_annotation | cut -f1 | uniq > $BAM_trim_dir/"chr_in_snpfile.txt"
rm -f $BAM_trim_dir/"snps_overlapping_annotation_resorted.bed"
while read line
do
	chr=`echo -e "$line" | cut -f1`
	in_snpfile=`awk -v c=$chr '$1==c' $BAM_trim_dir/"chr_in_snpfile.txt" | wc -l`
	if [ $in_snpfile != "0" ]; then # if more than 0 chr included
		awk -v OFS="\t" -v c=$chr '$1==c' $snps_overlapping_annotation >> $BAM_trim_dir/"snps_overlapping_annotation_resorted.bed"
	fi
done < $BAM_trim_dir/"chr_order_in_bam_header.txt"

intersectBed -abam $bam -b $BAM_trim_dir/"snps_overlapping_annotation_resorted.bed" -split -sorted -g $BAM_trim_dir/"chr_order_in_bam_header.txt" > $BAM_trim_dir"/overlap_reads_vs_snps.bam"

samtools view $BAM_trim_dir"/overlap_reads_vs_snps.bam" | perl -nae 'if ($F[1] & 0x10) {$strand = "minus" } else {$strand = "plus"}; $_ =~ s/$F[0]/${F[0]}_${strand}_${F[2]}_${F[3]}/; print $_;' | LANG=C sort -T $BAM_trim_dir > $BAM_trim_dir"/overlap_reads_vs_snps_for_join.sam"

echo -e "...create bedfile for reads overlapping snps."
bedtools bamtobed -bed12 -splitD -i $BAM_trim_dir"/overlap_reads_vs_snps.bam" | awk -v OFS="\t" '{gsub(/\/[12]/,"",$4); print $0}' > $BAM_trim_dir"/overlap_reads_vs_snps.bed"

echo -e "...intersect reads overlapping SNPs with SNPs overlapping annotation (for bed files)."
intersectBed -split -sorted -g $BAM_trim_dir/"chr_order_in_bam_header.txt" -a  $BAM_trim_dir"/overlap_reads_vs_snps.bed" -b $BAM_trim_dir/"snps_overlapping_annotation_resorted.bed" -wb -wa | awk -v OFS="\t" '{if($6=="+"){strand="plus";}else{strand="minus";} print $4"_"strand"_"$1"_"$2+1,$1,$2,$3,$6,$14,$16}' | LANG=C sort -k1 -T $BAM_trim_dir > $BAM_trim_dir"/overlap_data_reads_and_snps.txt"

echo -e "...divide reads overlapping SNPs in annotation into reads covering one SNP and reads covering multiple SNPs."
awk -v OFS=" " -v mult=$BAM_trim_dir"/mult_reads_single_line.txt" -v uniq=$BAM_trim_dir"/read_name.uniq" '{if($1==prev){snp=snp " " $7; pos=pos " " $6; info=$2" "$3" "$4" "$5; count++;} else { if (prev){ if(count>1){ print prev,count,info,snp,pos > mult; } else{ print prev > uniq;}} prev=$1; snp=$7; pos=$6; count=1;}} END {if (prev){ if(count>1){ print prev,count,info,snp,pos > mult;} else{ print prev > uniq;}}}' $BAM_trim_dir"/overlap_data_reads_and_snps.txt"

echo -e "...create uniq.sam file by joining."
samtools view -H $bam > $BAM_trim_dir"/uniq.sam"
LANG=C join -1 1 -2 1 $BAM_trim_dir"/overlap_reads_vs_snps_for_join.sam" $BAM_trim_dir"/read_name.uniq"  | awk -v OFS="\t" '{$1=$1};{print$0}' >> $BAM_trim_dir"/uniq.sam"

echo -e "...join file containing information about the multiple SNPs a read overlaps with the original sam file."
cut -f1-11 $BAM_trim_dir"/overlap_reads_vs_snps_for_join.sam" | LANG=C join -1 1 -2 1 - $BAM_trim_dir"/mult_reads_single_line.txt" | awk -v OFS="\t" '{$1=$1};{print$0}' > $BAM_trim_dir"/mult_reads_single_line_with_sam_data.txt"

echo -e "...create mult.sam using the perl script bamtrim.pl."
samtools view -H $bam > $BAM_trim_dir"/mult.sam"
cat "${BAM_trim_dir}/mult_reads_single_line_with_sam_data.txt" | perl $scripts_dir"/bamtrim.pl" >> $BAM_trim_dir"/mult.sam"

echo -e "...convert uniq.sam/mult.sam into .bam files."
samtools view -Sb $BAM_trim_dir"/uniq.sam" > $BAM_trim_dir"/uniq.bam"
samtools view -Sb $BAM_trim_dir"/mult.sam" > $BAM_trim_dir"/mult.bam"

echo -e "...create trimmed.bam by merging mult.bam and uniq.bam."
samtools merge -f $BAM_trim_dir"/trimmed.bam" $BAM_trim_dir"/mult.bam" $BAM_trim_dir"/uniq.bam"

echo -e "...sort trimmed.bam."
samtools sort $BAM_trim_dir"/trimmed.bam" -o $debug_dir"/trimmed_s.bam"
