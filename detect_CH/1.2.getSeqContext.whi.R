###
library(data.table)
## Seq context
library("BSgenome.Hsapiens.UCSC.hg19")
###

### command line input:
commandLineArgs <- commandArgs(TRUE) # 1.outdir 2./path/wl.csv.gz 3./path/dnmt.csv.gz 4./path/manul.csv.gz 5./path/outprefix 6./path/all_vars.csv.gz 7.batch_name 8.out_allvarfile


#############
outDir <- commandLineArgs[1]

wl_csv <- commandLineArgs[2]

dnmt_csv <- commandLineArgs[3]

manrev_csv <- commandLineArgs[4]

out_prefix <- commandLineArgs[5]

all_vars <- commandLineArgs[6] # 

batch_name <- commandLineArgs[7] # "whimips_longitudinal_20230412"

# out_allvarfile <- paste0(out_prefix,".",batch_name,".tsv")
out_allvarfile <- commandLineArgs[8]
  ##### 
wl.chip <- fread(wl_csv, header=T)
wl.chip$chr_pos_ref_alt <- paste(wl.chip$Otherinfo4,wl.chip$Otherinfo5, 
                                 wl.chip$Otherinfo7, wl.chip$Otherinfo8, 
                                 sep = "_")

dnmt_notwl <- fread(dnmt_csv, header=T)
dnmt_notwl$chr_pos_ref_alt <- paste(dnmt_notwl$Otherinfo4,dnmt_notwl$Otherinfo5, 
                                    dnmt_notwl$Otherinfo7, dnmt_notwl$Otherinfo8, 
                                 sep = "_")

man_review <- fread(manrev_csv, header=T)
man_review$chr_pos_ref_alt <- paste(man_review$Otherinfo4,man_review$Otherinfo5, 
                                    man_review$Otherinfo7, man_review$Otherinfo8, 
                                    sep = "_")
### 
all_var <- fread(all_vars, header=T)

all_var$chr_pos_ref_alt <- paste(all_var$Otherinfo4, all_var$Otherinfo5, 
                                 all_var$Otherinfo7, all_var$Otherinfo8, 
                                    sep = "_")
all_var$AF_raw <- as.numeric(all_var$AF_raw)
# only cosmic 96 (no cosmic70)
all_var[which(all_var$AF_raw>0.001),c(1,2,11,95,96:101,13,17)]

all_var <- subset(all_var, all_var$ExonicFunc.refGene!="synonymous SNV" & 
                    (all_var$AF_raw<0.001 | is.na(all_var$AF_raw)) )
                  
## Add seq context

chipVars <- all_var[,c(83,84,86,87, 101)]
names(chipVars) <- c("CHROM", "POS", "REF", 
                     "ALT", "chr_pos_ref_alt")

chipVars$CHROM <- paste0("chr",chipVars$CHROM)
chipVars$context=""



getflank <- function(position, alleles="[N/N]", chr="chr12", offset=10) {
  leftflank  <- getSeq(Hsapiens,chr,position-offset,position-1)
  rightflank <- getSeq(Hsapiens,chr,position+1,position+offset)
  paste(leftflank,alleles,rightflank,sep="")
}

for(i in 1:dim(chipVars)[1]){
  chipVars$context[i]= getflank(chipVars$POS[i], 
                                alleles=paste("[",chipVars$REF[i],"/",chipVars$ALT[i],"]",sep = ""),
                                chr=chipVars$CHROM[i]
  )
  if(i%%50==0){print(i)}
}

## Merge Sequence context to CHIP variant file
all_var <- merge(all_var, chipVars, by="chr_pos_ref_alt")

names(all_var)

my_all_var <- all_var[, c(1:3, 103:105, 8, 10, 11, 
                          93:95, 12, 17, 96:101,106)]

fwrite(my_all_var, out_allvarfile,
       row.names = F, col.names = T, 
       sep="\t", quote = F)

#########################
