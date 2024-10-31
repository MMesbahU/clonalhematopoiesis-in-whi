#!/bin/bash


source scripts/useuse

use VCFtools

use Bcftools
###################
out_path=${1} 

in_vcf=${2} 

varFile=${3} 

n_start=${4} # 1

n_end=${5} # 500

outFile="${out_path}/sample_${n_start}_${n_end}.vars_filt_min1_vaf1p0pct_5altrd.tsv"
#######################

echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > ${outFile}

while read samplename; do paste <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile}) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> ${outFile}; done < <(bcftools query -l ${in_vcf} | awk -v N_Start=${n_start} -v N_End=${n_end} 'NR>=N_Start && NR<=N_End')


gzip -f ${outFile}


