
## use script: "myWHI_20221123.script.sh"
## CHIP call: vaf>=0.1%; min forward/reverse Alt read>=2; min AD>=5; min DP>=50

###
library(data.table)
###
batch_number <- 20230711
###
setwd(paste0("whimips_longitudinal_",
             batch_number,
             "/per_sample_ch_var/"))
###

whi_chip <- fread("vaf001ad5dp50.vars_filt_min1_vaf1p0pct_5altrd.tsv.gz")

my_all_var <- fread(paste0("../all_putative_CHIP.whimips_longitudinal_",batch_number,".tsv"), 
                           header = T)

annot.whi_chip <- merge(whi_chip, 
                        my_all_var[,c(1,2,7:15,21)], 
				                by.x="CHROM_POS_REF_ALT", 
								        by.y="chr_pos_ref_alt")

												  
names(annot.whi_chip)

annot.whi_chip <- annot.whi_chip[,c(2,1,3:6,11,8:10,12:21,7)]

names(annot.whi_chip) <- c(names(annot.whi_chip)[1:17], 
                           "gnomAD_AF", "whitelist", 
                           "Seq_Context", "VEP_Annot")

# All putative CHIP variants
fwrite(annot.whi_chip, paste0("WHI_",batch_number,".chip_variants.vaf001_DP50_AD5_FR2.all_putative_CHIP.csv.gz"), 
       row.names = F, col.names = T, quote = T, sep = ",")

# Variants in White list Variants
fwrite(annot.whi_chip[annot.whi_chip$whitelist=="TRUE", ], 
          paste0("WHI_",batch_number,".chip_variants.vaf001_DP50_AD5_FR2.varOI_wl.csv.gz"), 
          row.names = F, col.names = T, quote = T, sep = ",")

# Other Variants not in (White list, DNMT3A variants + Manual Review Variants)
fwrite(annot.whi_chip[annot.whi_chip$whitelist=="FALSE",], 
  paste0("WHI_",batch_number,".chip_variants.vaf001_DP50_AD5_FR2.varOI_nonWL_vars.csv.gz"), 
  row.names = F, col.names = T, quote = T, sep = ",")


### exclude black listed variants provided by Sidd
var_blacklisted <- fread("whi/blacklist.csv", 
                         header=F)

# wl + no blacklist Variants
annot.whi_chip_wl_noBlacklist <- subset(annot.whi_chip, 
                                        !(annot.whi_chip$CHROM_POS_REF_ALT %in% var_blacklisted$V1) & 
                                          annot.whi_chip$whitelist=="TRUE")

fwrite(annot.whi_chip_wl_noBlacklist, 
          paste0("WHI_", batch_number,".chip_variants.vaf001_DP50_AD5_FR2.varOI_wl_noBlacklist.csv.gz"), 
          row.names = F, col.names = T, quote = T, sep = ",")

## all w/o blacklist
fwrite(annot.whi_chip[!(annot.whi_chip$CHROM_POS_REF_ALT %in% var_blacklisted$V1),], 
       paste0("WHI_",batch_number,".chip_variants.vaf001_DP50_AD5_FR2.all_putative_CHIP_noBlacklist.csv.gz"), 
       row.names = F, col.names = T, quote = T, sep = ",")


###### not-wl filter
head(sort(table(annot.whi_chip$CHROM_POS_REF_ALT
                [annot.whi_chip$VAF>=0.05 & 
                    annot.whi_chip$whitelist=="FALSE"]), 
          decreasing = T),20)

## only keep variants with "haematopoietic_and_lymphoid_tissue"
annot.whi_chip_wl_pls <- subset(annot.whi_chip, (annot.whi_chip$whitelist=="TRUE" | 
                                                   (annot.whi_chip$whitelist=="FALSE" & 
                                                      grepl(pattern = "haematopoietic_and_lymphoid_tissue",x = annot.whi_chip$cosmic96_coding) ) ) &
                                !(annot.whi_chip$CHROM_POS_REF_ALT %in% var_blacklisted$V1))

head(sort(table(annot.whi_chip_wl_pls$CHROM_POS_REF_ALT
                [annot.whi_chip_wl_pls$VAF>=0.05 & 
                    annot.whi_chip_wl_pls$whitelist=="FALSE"]), 
          decreasing = T),20)

fwrite(annot.whi_chip_wl_pls, 
       paste0("WHI_",batch_number,".chip_variants.vaf001_DP50_AD5_FR2.varOI_wlpls_noBlacklist.csv.gz"), 
       row.names = F, col.names = T, quote = T, sep = ",")
