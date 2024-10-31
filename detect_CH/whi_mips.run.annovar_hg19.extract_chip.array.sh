#!/bin/bash

#########################
## Mar 28, 2023
## qsub -wd whimips_longitudinal_20230313/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 1 -binding linear:1 -N vcf_annot_whi_20230313 detect_CH/whi_mips.run.annovar_hg19.extract_chip.array.sh annovar driver_genes.rda whimips_longitudinal_20230313/infoOnly.sorted.vars_filt_min1_vaf1pct_5altrd.4.vcf.gz whimips_longitudinal_20230313 whimips_longitudinal_20230313 "annot." "annot.infoOnly.sorted." ".4.hg19_multianno.txt.gz" detect_CH/run.CHIP_extract.update_29March2022.R
#########################

#########################
source /scripts/useuse

use Tabix

use Bcftools

use R-4.0

# use BEDTools
# use .perl-5.28.0
#######################

path_annovar=${1} # 

driver_genes=${2} # chip genes, rules

in_vcf_file=${3}

annot_out_dir=${4}

chip_out_dir=${5}

Prefix=${6}

suffix2move=${7} # annot.infoOnly.sorted.

prefix_2_remove=${8} # .4.hg19_multianno.txt.gz

r_script=${9}
##
# in_vcf_file=`awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${in_vcf_file_list} `

out_prefix=${annot_out_dir}/${Prefix}$(basename ${in_vcf_file} ".vcf.gz" )

############################################### Annovar ################################################

## Run annovar
${path_annovar}/table_annovar.pl ${in_vcf_file} ${path_annovar}/humandb/ \
	-buildver hg19 \
	--out ${out_prefix} \
	-remove \
	-protocol refGene,cosmic96_coding,gnomad211_exome,dbnsfp41a \
	-operation g,f,f,f \
	-nastring . \
	-vcfinput \
	-polish

# 
rm ${out_prefix}.avinput

bgzip ${out_prefix}.hg19_multianno.txt

bgzip ${out_prefix}.hg19_multianno.vcf

tabix -p vcf ${out_prefix}.hg19_multianno.vcf.gz

########################################################################################################
####### Run Rscript to extract CHIPs
        ## get CHIP
Rscript ${r_script} ${annot_out_dir} ${driver_genes} ${out_prefix}.hg19_multianno.txt.gz ${suffix2move} ${prefix_2_remove}

########################################################################################################

