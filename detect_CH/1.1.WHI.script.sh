
## 
## New Batch 
	## VCF file
use Google-Cloud-SDK

## Download vcf file from TERRA

## 2024-03-26

use Google-Cloud-SDK

BATCH="20240322"

mkdir -p whi/whimips_longitudinal_${BATCH}

medpop_outDir=whimips_longitudinal_${BATCH}
####

gsutil -m cp gs://gcp_path/whimips_longitudinal_${BATCH}/vars_filt_min1_vaf1p0pct_5altrd.5.vcf.gz ${medpop_outDir}/

gsutil -m cp gs://gcp_path/whimips_longitudinal_${BATCH}/varscan_calls_min1sampvaf1pct_min5altrd_mind50.txt.gz ${medpop_outDir}/

# bam in this batch

gsutil ls gs://gcp_path/whimips_longitudinal_${BATCH}/bam/*.bam >> ${medpop_outDir}/bam_list.${BATCH}.list

# all bam
gsutil ls gs://fc-2800aa50-51bc-487c-b599-64f3341c6265/whimips_longitudinal_*/bam/*.bam > ${medpop_outDir}/bam_list.whimips_longitudinal_all_till_${BATCH}.list


## no duplicate ids
echo -e "entity:mar2024_whimips_sample_id\tbam\tbai" > ${medpop_outDir}/bam_list_4_terra.till_${BATCH}.tsv

paste -d '\t' <(cat ${medpop_outDir}/bam_list.whimips_longitudinal_all_till_${BATCH}.list | sed 's:/:\t:g' | cut -f6 | sed 's:.bam::g') <(cat ${medpop_outDir}/bam_list.whimips_longitudinal_all_till_${BATCH}.list) <(awk '{print $1".bai"}' ${medpop_outDir}/bam_list.whimips_longitudinal_all_till_${BATCH}.list ) | awk '!seen[$1]++' >> ${medpop_outDir}/bam_list_4_terra.till_${BATCH}.tsv

##
cat bam_list_4_terra.till_20240322.tsv | grep -vE 'posctrl|opentrons_test_gDNA' | wc
# 14219
wc bam_list_4_terra.till_20240322.tsv    
# 14276



## modify VCF header
###### v.v.i to change:
### vi vcf to mannually change
## changed VCF info for "IS_POSCTRL_SITE" and "IS_IN_GNOMAD" from Number=. to Number=0
####

## add var column
batch_number=20240322
inFile='vars_filt_min1_vaf1p0pct_5altrd.5.vcf.gz'
zgrep '^#' whimips_longitudinal_${batch_number}/${inFile} >  whimips_longitudinal_${batch_number}/whi_${batch_number}.sorted.vars_filt_min1_vaf1pct_5altrd.vcf && paste -d'\t' <(zgrep -v '^#' whimips_longitudinal_${batch_number}/${inFile} | awk '{print $1"\t"$2"\t"$1"_"$2"_"$4"_"$5}' ) <(zgrep -v '^#' whimips_longitudinal_${batch_number}/${inFile} | cut -f4-) >> whimips_longitudinal_${batch_number}/whi_${batch_number}.sorted.vars_filt_min1_vaf1pct_5altrd.vcf &

###### v.v.i to change:
### vi vcf to mannually change
## changed VCF info for "IS_POSCTRL_SITE" and "IS_IN_GNOMAD" from Number=. to Number="A"

# gzip whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf
sort_input=whimips_longitudinal_${batch_number}/whi_${batch_number}.sorted.vars_filt_min1_vaf1pct_5altrd.vcf

gzip ${sort_input}

##### Samples
use Bcftools
# bcftools -v
# bcftools 1.16
# Using htslib 1.16
# Copyright (C) 2022 Genome Research Ltd.
bcftools query -l ${sort_input}.gz > $(dirname ${sort_input}.gz)/whi_${batch_number}_sample.list


# info only vcf
vcf_in=whimips_longitudinal_${batch_number}/whi_${batch_number}.sorted.vars_filt_min1_vaf1pct_5altrd.vcf

zgrep '^##' ${vcf_in} > $(dirname ${vcf_in})/infoOnly.$(basename ${vcf_in} ".gz") && echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> $(dirname ${vcf_in})/infoOnly.$(basename ${vcf_in} ".gz") &&  zgrep -v '^#' ${vcf_in} | cut -f1-5 | awk '{print $0"\t.\t.\t."}' >> $(dirname ${vcf_in})/infoOnly.$(basename ${vcf_in} ".gz") &

gzip $(dirname ${vcf_in})/infoOnly.$(basename ${vcf_in} ".gz")


#### Run annovar
## Copy annovar files
## get CHIP
## updated, current:
qsub -wd whimips_longitudinal_${batch_number}/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 1 -binding linear:1 -N vcf_annot_whi_${batch_number} detect_CH/whi_mips.run.annovar_hg19.extract_chip.array.sh annovar driver_genes.rda $(dirname ${vcf_in})/infoOnly.$(basename ${vcf_in} ".gz").gz $(dirname ${vcf_in}) $(dirname ${vcf_in}) "annot." "infoOnly.whi_${batch_number}.sorted." ".hg19_multianno.txt.gz" detect_CH/run.CHIP_extract.update_29March2022.R 

	## Get Sequence Context
	## in Mac or run Rscript 
# whi_20221123/extractCHIP.whi_20221123.R
# use R-4.1
# detect_CH/extractCHIP.whi_20221123.R
batch_number=20240322; qsub -wd whimips_longitudinal_${batch_number}/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 1 -binding linear:1 -N get_Seq_context.whi_${batch_number} detect_CH/1.2.run.getSeqContext.sh detect_CH/1.2.getSeqContext.whi.R whimips_longitudinal_${batch_number} whimips_longitudinal_${batch_number}/annot.vars_filt_min1_vaf1pct_5altrd.varsOI.wl.csv.gz whimips_longitudinal_${batch_number}/annot.vars_filt_min1_vaf1pct_5altrd.varsOI.DNMT3Amis.notwl.csv.gz whimips_longitudinal_${batch_number}/annot.vars_filt_min1_vaf1pct_5altrd.varsOI.manualreview.csv.gz whimips_longitudinal_${batch_number}/all_putative_CHIP whimips_longitudinal_${batch_number}/annot.vars_filt_min1_vaf1pct_5altrd.annovar.varsOI.allvariants.csv.gz whimips_longitudinal_${batch_number}

####### Per Sample Variants ###########
### Get Variant info for each Sample
batch_number=20240322
mkdir -p whimips_longitudinal_${batch_number}/per_sample_ch_var
mkdir -p whi/whimips_longitudinal_${batch_number}/tmpdir

for i in `seq 1 500 $(wc -l whimips_longitudinal_${batch_number}/whi_${batch_number}_sample.list | awk '{print $1}')`; do qsub -R y -wd whimips_longitudinal_${batch_number}/tmpdir -N var_per_sam_${i}_$(( ${i} + 499 )) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 detect_CH/1.3.extract_CH_info.qsub.sh whimips_longitudinal_${batch_number}/per_sample_ch_var whimips_longitudinal_${batch_number}/whi_${batch_number}.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz whimips_longitudinal_${batch_number}/all_putative_CHIP.whimips_longitudinal_${batch_number}.tsv ${i} $(( ${i} + 499 )); done


####### Filter variants 
## 
######### VAF>=0.1%; minAD>=5; min FR&RR >=2; DP>=50
# echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/vaf001.whi_20221123.all_variants.vars_filt_dp50_min1_vaf1pct_5altrd.tsv

# for files in $(ls -lv whi_20221123/per_sample_ch_var/*.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv.gz | awk '{print $NF}'); do zcat ${files} | awk '(NR>1 && $8>=50 && $9>=0.001){split($10, DP4, ","); {if(DP4[3]>1 && DP4[4]>1 && (DP4[3]+DP4[4])>=5 ) print $0 }}' >> whi_20221123/vaf001.whi_20221123.all_variants.vars_filt_dp50_min1_vaf1pct_5altrd.tsv; done &

batch_number=20240322

outDir=whimips_longitudinal_${batch_number}/per_sample_ch_var

infileSuffix=".vars_filt_min1_vaf1p0pct_5altrd.tsv.gz"

outFilePrefix="vaf001ad5dp50"

echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4\tREF_FR\tREF_RR\tALT_FR\tALT_RR\tAD_REF\tAD_ALT" > ${outDir}/${outFilePrefix}.vars_filt_min1_vaf1p0pct_5altrd.tsv

for files in $(ls -lv ${outDir}/*${infileSuffix} | awk '{print $NF}'); do zcat ${files} | awk '(NR>1 && $8>=50 && $9>=0.001){split($10, DP4, ","); {if(DP4[3]>1 && DP4[4]>1 && (DP4[3]+DP4[4])>=5 ) print $0"\t"DP4[1]"\t"DP4[2]"\t"DP4[3]"\t"DP4[4]"\t"DP4[1]+DP4[2]"\t"DP4[3]+DP4[4]}}' >> ${outDir}/${outFilePrefix}.vars_filt_min1_vaf1p0pct_5altrd.tsv; done &

gzip -f ${outDir}/${outFilePrefix}.vars_filt_min1_vaf1p0pct_5altrd.tsv &

## all other variants with DP<50 or DP==0
outFilePrefix_v2="dp_less50"

echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4\tREF_FR\tREF_RR\tALT_FR\tALT_RR\tAD_REF\tAD_ALT" > ${outDir}/${outFilePrefix_v2}.vars_filt_min1_vaf1p0pct_5altrd.tsv

for files in $(ls -lv ${outDir}/*${infileSuffix} | awk '{print $NF}'); do zcat ${files} | awk '(NR>1 && $8<50 || $8==0 || $9<0.001){split($10, DP4, ","); {if(DP4[3]<=1 || DP4[4]<=1 || (DP4[3]+DP4[4])<5 ) print $0"\t"DP4[1]"\t"DP4[2]"\t"DP4[3]"\t"DP4[4]"\t"DP4[1]+DP4[2]"\t"DP4[3]+DP4[4]}}' >> ${outDir}/${outFilePrefix_v2}.vars_filt_min1_vaf1p0pct_5altrd.tsv; done &

gzip -f ${outDir}/${outFilePrefix_v2}.vars_filt_min1_vaf1p0pct_5altrd.tsv &

## variants with DP>=50 & VAF<0.001
outFilePrefix_v3="vafGrt001dp50"
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > ${outDir}/${outFilePrefix_v3}.vars_filt_min1_vaf1p0pct_5altrd.tsv

for files in $(ls -lv ${outDir}/*${infileSuffix} | awk '{print $NF}'); do zcat ${files} | awk '(NR>1 && $8>=50 && $9<0.001) {print $0 }' >> ${outDir}/${outFilePrefix_v3}.vars_filt_min1_vaf1p0pct_5altrd.tsv; done &

gzip -f ${outDir}/${outFilePrefix_v3}.vars_filt_min1_vaf1p0pct_5altrd.tsv &


## Add annovar annotations:
### next step in 
# R script "detect_CH/annot.annovar.R"
######### 

########
