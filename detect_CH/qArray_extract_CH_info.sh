#!/bin/bash


## all jobs
## for i in `seq 1 500 $(wc -l whi_20230412_sample.list | awk '{print $1}')`; do qsub -R y -wd whimips_longitudinal_20230412/tmpdir -N var_per_sam_${i}_$(( ${i} + 499 )) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 script/detect_CH/qArray_extract_CH_info.sh ${i} $(( ${i} + 499 )); done


source /scripts/useuse

use VCFtools

use Bcftools


out_path="whimips_longitudinal_20230412/per_sample_ch_var"

in_vcf="whimips_longitudinal_20230412/whi_20230412.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz"

varFile="whimips_longitudinal_20230412/all_putative_CHIP.whimips_longitudinal_20230412.tsv"

n_start=${1} # 1

n_end=${2} # 500

outFile="${out_path}/sample_${n_start}_${n_end}.vars_filt_min1_vaf1p0pct_5altrd.tsv"
##########
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > ${outFile}

while read samplename 
do 
    paste <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile} ) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf ${in_vcf} --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' ${varFile}) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> ${outFile} 

done < <(bcftools query -l ${in_vcf} | awk -v N_Start=${n_start} -v N_End=${n_end} 'NR>=N_Start && NR<=N_End')



