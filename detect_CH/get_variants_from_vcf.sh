#!/bin/bash


source scripts/useuse


use VCFtools

use Bcftools


### 
out_file=${1}

in_vcf_file=${2}
### Sample 1 500
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/sample_1_500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> whi_20221123/sample_1_500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=1 && NR<=500') &

### 
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > ${out_file}

## 
while read samplename; do paste <(vcftools --gzvcf ${in_vcf_file} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> whi_20221123/sample_1_500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=1 && NR<=500')

