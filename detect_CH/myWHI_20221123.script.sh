
## Sept 10, 2022
## New Batch 
	## VCF file
use Google-Cloud-SDK

gsutil cp gs://filepath/vars_filt_min1_vaf1pct_5altrd.4.vcf.gz .

## add var
zgrep '^#' whi_20221123/vars_filt_min1_vaf1pct_5altrd.4.vcf.gz >  whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf && paste -d'\t' <(zgrep -v '^#' whi_20221123/vars_filt_min1_vaf1pct_5altrd.4.vcf.gz | awk '{print $1"\t"$2"\t"$1"_"$2"_"$4"_"$5}' ) <(zgrep -v '^#' whi_20221123/vars_filt_min1_vaf1pct_5altrd.4.vcf.gz | cut -f4-) >> whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf &

## changed VCF info for "IS_POSCTRL_SITE" and "IS_IN_GNOMAD" from Number=. to Number="A"

gzip whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf

##### Samples
use Bcftools

bcftools query -l whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz > whi_20221123/whi_20221123_samples.list

# info only vcf
zgrep '^##' whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz > whi_20221123/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf && echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> whi_20221123/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf &&  zgrep -v '^#' whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | cut -f1-5 | awk '{print $0"\t.\t.\t."}' >> whi_20221123/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf &

gzip whi_20221123/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf


#### Run annovar
qsub -wd whi_20221123 -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 1 -binding linear:1 -N vcf_annot_whi_20221123 whi_20220823.run.annovar_hg19.extract_chip.array.sh annovar whi_20221123/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf.gz whi_20221123 whi/whi_20221123 "annot" ".hg19_multianno.txt.gz"
	## get CHIP
Rscript run.CHIP_extract.update_29March2022.R whi_20221123 whi_20221123/annotinfoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.hg19_multianno.txt.gz annotinfoOnly.sorted. .4.hg19_multianno.txt.gz
	## in Mac run rscript 
whi_20221123/extractCHIP.whi_20221123.R

### Get Variant info for each Sample
# echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv

use VCFtools
use Bcftools



# Sample 1 500
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/sample_1_500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> whi_20221123/sample_1_500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=1 && NR<=500') &

# Sample 501:1000
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/sample_501_1000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >> whi_20221123/sample_501_1000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=501 && NR<=1000') &

# Sample 1001_1500
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" > whi_20221123/sample_1001_1500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}' whi_20221123 all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >>    whi_20221123/sample_1001_1500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=1001 && NR<=1500') &

# Sample 1501_2000
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" >    whi_20221123/sample_1501_2000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >>    whi_20221123/sample_1501_2000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=1501 && NR<=2000') &

# Sample 2001_2500
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" >    whi_20221123/sample_2001_2500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >>    whi_20221123/sample_2001_2500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=2001 && NR<=2500') &

# Sample 2501_3000
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" >    whi_20221123/sample_2501_3000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >>    whi_20221123/sample_2501_3000.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=2501 && NR<=3000') &

# Sample 3001_3500
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" >    whi_20221123/sample_3001_3500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; while read samplename; do paste <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --get-INFO ANN | cut -f1-5 | awk -v mysamp=${samplename} 'NR>1{print mysamp"\t"$1"_"$2"_"$3"_"$4"\t"$0}') <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info AF | awk 'NR>1{print $3}' ) <(vcftools --gzvcf    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz --indv ${samplename} --snps <(awk -F'\t' 'NR>1{print $1}'    whi_20221123/all_putative_CHIP.whi_20221123.tsv) -c --extract-FORMAT-info DP4 | awk 'NR>1{print $3}')  >>    whi_20221123/sample_3001_3500.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv; done < <(bcftools query -l    whi_20221123/whi_20221123.sorted.vars_filt_min1_vaf1pct_5altrd.vcf.gz | awk 'NR>=3001 && NR<=3500') &



## zcat    whi_20221123/per_sample_ch_var/whi_plate45_WHI_CHIP_plate45_A10.whi_20221123.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv.gz | awk '$8>=50 && $9>=0.001{split($10, DP4, ","); {if(DP4[3]>1 && DP4[4]>1 && (DP4[3]+DP4[4])>=5 ) print $0 }}' | less -S

######### VAF>=0.1%; minAD>=5; min FR&RR >=2; DP>=50
echo -e "Sampleid\tCHROM_POS_REF_ALT\tCHROM\tPOS\tREF\tALT\tVEP_Annot\tDP\tVAF\tDP4" >    whi_20221123/vaf001.whi_20221123.all_variants.vars_filt_dp50_min1_vaf1pct_5altrd.tsv

for files in $(ls -lv    whi_20221123/per_sample_ch_var/*.all_variants.vars_filt_min1_vaf1pct_5altrd.tsv.gz | awk '{print $NF}'); do zcat ${files} | awk '(NR>1 && $8>=50 && $9>=0.001){split($10, DP4, ","); {if(DP4[3]>1 && DP4[4]>1 && (DP4[3]+DP4[4])>=5 ) print $0 }}' >> whi_20221123/vaf001.whi_20221123.all_variants.vars_filt_dp50_min1_vaf1pct_5altrd.tsv; done & 

### next step in R script "whi_20221123/extractCHIP.whi_20221123.R"
#####

