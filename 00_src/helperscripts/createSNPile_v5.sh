#!/bin/bash
## Script to generate a SNPfile .bed format as required by Allelome.PRO for v5
## Script was updated to fit the mgp.v5.merged.snps_all.dbSNP142.vcf
## Update l.19: 1_q[6] is for v3 vcf; the FI filter if variant passes filter = 1 -> in v5 vcf this is at POS 1_q[14]
## Tim Hasenbeim, 01.06.21
grep -v "^##" $1 | head -n1 | awk '{for(i=10;i<=NF;i++){print i-9": "$i}}'

echo "choose strain 1"
read column1
echo "choose strain 2"
read column2

((col1_in_file=column1+9))
((col2_in_file=column2+9))

# first check which element contains the information about whether the call passed the filter or not (FI)
grep -v "^#" $1 | awk -v OFS="\t" -v c1=$col1_in_file -v c2=$col2_in_file '{n=split($9,f,":"); if(f[n]=="FI"){fi=n;} else{ for(i in f){if(f[i]=="FI"){fi=i; break;}}} split(substr($c1,1,3),s1_i,"/"); split(substr($c2,1,3),s2_i,"/"); if(s1_i[1]==s1_i[2] && s2_i[1]==s2_i[2] && s1_i[1]!=s2_i[1]){ \
split($c1,s1_q,":"); split($c2,s2_q,":"); \
if(s1_q[14]==1 && s2_q[14]==1) { split($5,var,",");  \
if(s1_i[1]==0){s1_var=$4;}else{s1_var=var[s1_i[1]];} if(s2_i[1]==0){s2_var=$4;}else{s2_var=var[s2_i[1]];} \
print "chr"$1,$2,$2+1,s1_var""s2_var,"0","."}}}' | sort -k1,1 -k2,2n > SNPfile.bed
