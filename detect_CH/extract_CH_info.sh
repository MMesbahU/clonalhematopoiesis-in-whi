#!/bin/bash



source scripts/useuse

use VCFtools

# use Bcftools
echo $date

row_start=${1}
row_end=${2}
outDir=${3} # 
input_vcf=${4} 
ch_vars=${5} 
sample_list=${6} 
## total sample: 3313
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > ${outDir}/sample_${row_start}_${row_end}.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv

while read samplename; do echo ${samplename}; paste <(vcftools --gzvcf ${outDir}/${input_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${ch_vars} ) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf ${outDir}/${input_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${ch_vars} ) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${outDir}/${input_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${ch_vars} ) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${outDir}/${input_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${ch_vars} ) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}') |  awk '$9>0' >> ${outDir}/sample_${row_start}_${row_end}.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(awk -v ROW_start=${ROW_start} -v ROW_end=${ROW_end} 'NR>=ROW_start && NR<=ROW_end' ${sample_list})

## 
gzip ${outDir}/sample_${row_start}_${row_end}.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv




