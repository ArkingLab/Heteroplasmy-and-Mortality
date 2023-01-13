## For generating het / homo count, per sample and per pos file
install.packages("tidyverse")
install.packages("data.table")

## Load
library(tidyverse)
library(data.table)

## Load in files from MitoHPC output
#counts = read.table('mutect2.mutect2.05.all.0516.vcf', skip = 52)
counts = read.table('mutect2.mutect2.05.vcf', skip = 60)
colnames(counts) = c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE')
counts$SAMPLE = as.numeric(substr(counts$SAMPLE,1,7))

## Filter for variant level FILTER
counts = counts %>% filter(!str_detect(FILTER,"base|strand|slippage|weak|germline|position"))

## Filter for variant level INFO
counts = counts[!grepl('INDEL',counts$INFO),]
counts = counts[!grepl('HP',counts$INFO),]
counts1 = counts[!grepl('NUMT',counts$INFO),]
counts2 = counts[grepl('NUMT',counts$INFO),]
counts1$filter_numt = 'notnumt'
counts2$filter_numt = 'numt'
counts = rbind(counts1,counts2)
rm(counts1)
rm(counts2)

## Get info from INFO column
test = unlist(strsplit(counts$INFO, split=';'))
AF = test[grepl('AF',test)]
GT = test[grepl('GT',test)]
DP = test[grepl('DP',test)]
HG = test[grepl('HG',test)]
NONSYN = test[grepl('NONSYN',test)]
STOP = test[grepl('STOP',test)]
AP = test[grepl('AP',test)]
APS = test[grepl('APS',test)]
counts$AF = as.numeric(substr(AF,4,15))
counts$Genotype = substr(GT,4,15)
counts$Read_depth = as.numeric(substr(DP,4,15))
counts$Haplogroup = NA

## Filter for Read_counts < 300
counts = counts[counts$Read_depth >= 300,]

## Get info from INFO column
test = unlist(strsplit(counts$INFO, split=';'))
HG = test[grepl('HG',test)]
NONSYN = test[grepl('NONSYN',test)]
STOP = test[grepl('STOP',test)]
AP = test[grepl('AP',test)]
APS = test[grepl('APS',test)]

## haplogroup
haplogroup = read.table('MitoHPC/mutect2.haplogroup1.tab',header=TRUE)
counts = merge(counts,haplogroup,by='SAMPLE')

## Get genes and complex information
countstest = counts[!grepl('DLOOP',counts$INFO),]
countstest = countstest[!grepl('TRN',countstest$INFO),]
countstest = countstest[!grepl('RNR',countstest$INFO),]
countstest = countstest[!grepl('CDS',countstest$INFO),]
countstest$GENE = NA
countstest$COMPLEX = NA
counts$GENE = NA
counts$COMPLEX = NA

countstrn = counts[grepl('TRN',counts$INFO),]
trnlist = unlist(strsplit(countstrn$INFO, split=';'))
trngene = trnlist[grepl('TRN',trnlist)]
countstrn$GENE = substr(trngene,5,15)
countstrn$COMPLEX = 'TRNA'

countscds = counts[grepl('CDS',counts$INFO),]
cdslist = unlist(strsplit(countscds$INFO, split=';'))
cdsgene = cdslist[grepl('CDS',cdslist)]
countscds$GENE = substr(cdsgene,5,15)
cdscomplex = cdslist[grepl('COMPLEX',cdslist)]
countscds$COMPLEX = substr(cdscomplex,9,15)

countsrnr = counts[grepl('RNR=',counts$INFO),]
rnrlist = unlist(strsplit(countsrnr$INFO, split=';'))
rnrgene = rnrlist[grepl('RNR=',rnrlist)]
countsrnr$GENE = substr(rnrgene,5,15)
countsrnr$COMPLEX = 'RRNA'

countsdloop = counts[grepl('DLOOP',counts$INFO),]
countsdloop$GENE = 'DLOOP'
countsdloop$COMPLEX = 'DLOOP'

countsnew = rbind(countstrn, countscds)
countsnew = rbind(countsnew, countsrnr)
countsnew = rbind(countsnew, countsdloop)
countsnew = rbind(countsnew, countstest)

counts = countsnew
rm(countsnew)
rm(countscds)
rm(countsdloop)
rm(countsrnr)
rm(countstest)
rm(countstrn)

## unique
counts$unique = paste(paste(counts$POS, counts$REF, sep = '_'), counts$ALT, sep = '_')

## new mito score
scores1 = fread('new_mito_score.tsv')
scores1only = scores1
scores1only$unique = paste(paste(scores1only$POS, scores1only$REF, sep = '_'), scores1only$ALT, sep = '_')
scores1only = scores1only[!duplicated(scores1only$unique),]
scores1only = scores1only[,c(1:3,5:6)]
countsscore = merge(counts, scores1only, by = c('POS','REF','ALT'), all.x = TRUE)
colnames(countsscore) = c(colnames(countsscore)[1:19],'mito_lc_consequence','mito_lc_score')
counts = countsscore
rm(countsscore)
rm(scores1)
rm(scores1only)

## Other scores
test = unlist(strsplit(counts$INFO, split=';'))
AP = test[grepl('AP',test)]
APS = test[grepl('APS',test)]
MMC = test[grepl('MMC',test)]
MLC = test[grepl('MLC',test)]
MCC = test[grepl('MCC',test)]

## Apogee score
AP_label = AP[seq(1,length(AP),by=2)]
AP_score = AP[seq(2,length(AP),by=2)]
counts$AP_label = NA
counts$AP_score = NA
countsap1 = counts[grepl('APS',counts$INFO),]
countsap2 = counts[!grepl('APS',counts$INFO),]
countsap1$AP_label = substr(AP_label,4,500)
countsap1$AP_score = as.numeric(substr(AP_score,5,20))
counts = rbind(countsap1,countsap2)
rm(countsap1)
rm(countsap2)

## Missense score
test = unlist(strsplit(counts$INFO, split=';'))
MMC = test[grepl('MMC',test)]
MMC_score = MMC[seq(2,length(MMC),by=3)]
MMC_class = MMC[seq(1,length(MMC),by=3)]
MMC_consq = MMC[seq(3,length(MMC),by=3)]
counts$MMC_score = NA
counts$MMC_class = NA
counts$MMC_consq = NA
countsmmc1 = counts[grepl('MMC',counts$INFO),]
countsmmc2 = counts[!grepl('MMC',counts$INFO),]
countsmmc1$MMC_score = as.numeric(substr(MMC_score,11,50))
countsmmc1$MMC_class = substr(MMC_class,11,50)
countsmmc1$MMC_consq = substr(MMC_consq,11,50)
counts = rbind(countsmmc1,countsmmc2)
rm(countsmmc1)
rm(countsmmc2)

## protein_gene_missense_constraint; missense_OEUF (MCC)
test = unlist(strsplit(counts$INFO, split=';'))
MCC = test[grepl('MCC',test)]
counts$MCC_score = NA
countsmcc1 = counts[grepl('MCC',counts$INFO),]
countsmcc2 = counts[!grepl('MCC',counts$INFO),]
countsmcc1$MCC_score = as.numeric(substr(MCC,5,20))
counts = rbind(countsmcc1,countsmcc2)
rm(countsmcc1)
rm(countsmcc2)

## Sep homo and het variants
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]

## Filter out those with two ALT allele that add up to 1 (biallelic het)
countshetbi = countshet %>% group_by(SAMPLE,POS) %>% slice_max(mito_lc_score, with_ties = FALSE) %>% ungroup()
countshet = countshetbi
rm(countshetbi)

## Counts of het and homo
uniquecountshet = sort(table(countshet$unique),decreasing = TRUE)
uniquecountshomo = sort(table(countshomo$unique),decreasing = TRUE)
uniquecountshetdf = as.data.frame(uniquecountshet)
uniquecountshomodf = as.data.frame(uniquecountshomo)
countshetwithfreq = merge(countshet, uniquecountshetdf, by.x = 'unique', by.y = 'Var1')
countshomowithfreq = merge(countshomo, uniquecountshomodf, by.x = 'unique', by.y = 'Var1')
countshet = countshetwithfreq
countshomo = countshomowithfreq
rm(countshetwithfreq)
rm(countshomowithfreq)
colnames(countshet) = c(colnames(countshet)[1:26], 'Freq_Het')
colnames(countshomo) = c(colnames(countshomo)[1:26], 'Freq_Homo')

## Merge freqs
countshetfreqonly = countshet[,c(1,27)]
countshetfreqonly = countshetfreqonly[!duplicated(countshetfreqonly$unique),]
countshomofreqonly = countshomo[,c(1,27)]
countshomofreqonly = countshomofreqonly[!duplicated(countshomofreqonly$unique),]
countshetwithhomo = merge(countshet, countshomofreqonly, by = 'unique', all.x = TRUE)
countshomowithhet = merge(countshomo, countshetfreqonly, by = 'unique', all.x = TRUE)
countshomowithhet = countshomowithhet[,c(1:26,28,27)]
countshet = countshetwithhomo
countshomo = countshomowithhet
rm(countshetwithhomo)
rm(countshomowithhet)
rm(countshetfreqonly)
rm(countshomofreqonly)

## Calculate Median and Max
counts = countshet
allunique = unique(counts$unique)
medianper = counts[counts$unique == allunique[1],]
medianper$Median_AF_Het = median(medianper$AF)
medianper$Max_AF_Het = max(medianper$AF)
newtotalmediantable = medianper
for (i in 2:length(allunique)) {
  medianper = counts[counts$unique == allunique[i],]
  medianper$Median_AF_Het = median(medianper$AF)
  medianper$Max_AF_Het = max(medianper$AF)
  newtotalmediantable = bind_rows(newtotalmediantable, medianper)
}
countshet = newtotalmediantable
rm(newtotalmediantable)
countshomo$Median_AF_Het = NA
countshomo$Max_AF_Het = NA
counts = rbind(countshet, countshomo)

## Transition and Transversion
counts$mutation = paste(counts$REF, counts$ALT, sep = 'to')
counts$Substitution = ifelse((counts$mutation == 'AtoG' | counts$mutation == 'GtoA' | counts$mutation == 'CtoT' | counts$mutation == 'TtoC'),'Transition','Transversion')
counts1 = counts[counts$Substitution == 'Transition',]
counts2 = counts[counts$Substitution == 'Transversion',]
counts2$Substitution = ifelse((counts2$mutation %in% c('AtoC','CtoA','AtoT','TtoA','GtoC','CtoG','GtoT','TtoG')),'Transversion',NA)
counts = rbind(counts1, counts2)

## Non-syn
counts$mutation_nonsynonymous = NA
counts$mutation_stop = NA
countsnons1 = counts[grepl('NONSYN',counts$INFO),]
countsnons2 = counts[!grepl('NONSYN',counts$INFO),]
countsnons1$mutation_nonsynonymous = 'NONSYN'
countsnons = rbind(countsnons1,countsnons2)
counts = countsnons

## STOP
countsstop1 = counts[grepl('STOP',counts$INFO),]
countsstop2 = counts[!grepl('STOP',counts$INFO),]
countsstop1$mutation_stop = 'STOP'
countsstop = rbind(countsstop1,countsstop2)
counts = countsstop

## Sep and id
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]
countshomo = transform(countshomo, identifier=as.numeric(factor(countshomo$unique)))
countshet = transform(countshet, identifier=as.numeric(factor(countshet$unique)))

## Download version (no genotype information)
countshomodownload = countshomo[,c(8,11:30,32:35)]
countshetdownload = countshet[,c(8,11:30,32:35)]

## Save files
write.table(countshomo, file='all_variants_homo_allfilter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshomodownload, file='all_variants_homo_download_allfilter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshet, file='all_variants_het_allfilter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshetdownload, file='all_variants_het_download_allfilter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)

## map file
het_id_unique_map = countshet[!duplicated(countshet$unique),c(1,35)]
homo_id_unique_map = countshomo[!duplicated(countshomo$unique),c(1,35)]
write.table(het_id_unique_map, file='het_id_unique_map.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(homo_id_unique_map, file='homo_id_unique_map.txt', sep="\t", row.names=FALSE, quote = FALSE)

## Read in files
countshet = fread('all_variants_het_allfilter_0621.txt')
persample = fread('per_sample_count_with_all_samplefilter_0621.txt')

## Per sample counts
allsample = unique(counts$SAMPLE)
usefullper = counts[counts$SAMPLE == allsample[1],]
usefullper$count_total = nrow(usefullper)
usefullperhet = usefullper[!usefullper$AF == 1,]
usefullperhomo = usefullper[usefullper$AF == 1,]
usefullper$count_het = nrow(usefullperhet)
usefullper$count_homo = nrow(usefullperhomo)
newtotaltable = usefullper[1,c(11,35:37)]

## for loop
for (i in 2:length(allsample)) {
  usefullper = counts[counts$SAMPLE == allsample[i],]
  usefullper$count_total = nrow(usefullper)
  usefullperhet = usefullper[!usefullper$AF == 1,]
  usefullperhomo = usefullper[usefullper$AF == 1,]
  usefullper$count_het = nrow(usefullperhet)
  usefullper$count_homo = nrow(usefullperhomo)
  newtotaltable = bind_rows(newtotaltable, usefullper[1,c(11,35:37)])
}

## Add in missed samples
everyone = fread('cvg.tab')
everyone$SAMPLE = as.numeric(substr(everyone$Run,1,7))
allsamplesmissed = everyone$SAMPLE[!everyone$SAMPLE %in% newtotaltable$SAMPLE]
missedsampletable = data.frame(SAMPLE = allsamplesmissed, count_total = 0, count_het = 0, count_homo = 0)
newtotaltable = rbind(newtotaltable, missedsampletable)

## Add in filter criteria for samples
sus1 = fread('mutect2.suspicious.tab')
sus1$SAMPLE = as.numeric(substr(sus1$V1,1,7))
sus11 = sus1[sus1$V2 == 'mean_cvg_less_500',]
sus12 = sus1[sus1$V2 == 'min_cvg_less_100',]
newtotaltable$mean_cvg_less_500 = ifelse((newtotaltable$SAMPLE %in% sus11$SAMPLE),'YES','NO')
newtotaltable$min_cvg_less_100 = ifelse((newtotaltable$SAMPLE %in% sus12$SAMPLE),'YES','NO')

sus2 = fread('mutect2.05.suspicious.tab1', header = FALSE)
sus2$SAMPLE = as.numeric(substr(sus2$V1,1,7))
sus21 = sus2[sus2$V2 == 'haplocheck_fail',]
sus22 = sus2[sus2$V2 == 'mismatch_HG',]
sus23 = sus2[sus2$V2 == 'multiple_NUMTs',]
newtotaltable$haplocheck_fail = ifelse((newtotaltable$SAMPLE %in% sus21$SAMPLE),'YES','NO')
newtotaltable$mismatch_HG = ifelse((newtotaltable$SAMPLE %in% sus22$SAMPLE),'YES','NO')
newtotaltable$multiple_NUMTs = ifelse((newtotaltable$SAMPLE %in% sus23$SAMPLE),'YES','NO')

sus3 = fread('exclude.list_May2022')
sus31 = sus3[sus3$reason == 'contaminated',]
sus32 = sus3[sus3$reason == 'low_CN',]
newtotaltable$excludelist_contaminated = ifelse((newtotaltable$SAMPLE %in% sus31$subject),'YES','NO')
newtotaltable$excludelist_lowCN = ifelse((newtotaltable$SAMPLE %in% sus32$subject),'YES','NO')
newtotaltable$hetp_count_greater5 = ifelse((newtotaltable$count_het > 5),'YES','NO')

cn = fread('WGS_CN.txt')
newtotaltable1 = merge(newtotaltable, cn, by.x = 'SAMPLE', by.y = 'subject', all.x = TRUE)
newtotaltable1[is.na(newtotaltable1)] = 0
newtotaltable = newtotaltable1
newtotaltable$lowCN = ifelse((newtotaltable$CN <= 40),'YES','NO')

write.table(newtotaltable, file='per_sample_count_with_all_filter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)

## Per position counts
## Filter for samples passed all filters only
filteredsample = newtotaltable %>%
  filter(mean_cvg_less_500=='NO',min_cvg_less_100=='NO',haplocheck_fail=='NO',mismatch_HG=='NO',multiple_NUMTs=='NO',
         excludelist_contaminated=='NO',excludelist_lowCN=='NO',hetp_count_greater5=='NO',lowCN=='NO')
countspos = counts[counts$SAMPLE %in% filteredsample$SAMPLE,]

allpos = unique(countspos$POS)
allunique = unique(countspos$unique)
usefullper = countspos[countspos$unique == allunique[1],]
usefullper$count_total = nrow(usefullper)
usefullperhet = usefullper[!usefullper$AF == 1,]
usefullperhomo = usefullper[usefullper$AF == 1,]
usefullper$count_het = nrow(usefullperhet)
usefullper$count_homo = nrow(usefullperhomo)
newpostable = usefullper[1,c(1:4,35:37,17:26,29:34)]

## for loop
for (i in 2:length(allunique)) {
  usefullper = countspos[countspos$unique == allunique[i],]
  usefullper$count_total = nrow(usefullper)
  usefullperhet = usefullper[!usefullper$AF == 1,]
  usefullperhomo = usefullper[usefullper$AF == 1,]
  usefullper$count_het = nrow(usefullperhet)
  usefullper$count_homo = nrow(usefullperhomo)
  newpostable = bind_rows(newpostable, usefullper[1,c(1:4,35:37,17:26,29:34)])
}
newpostable = newpostable[order(newpostable$POS),]
write.table(newpostable, file='per_pos_ref_alt_count_with_all_filter_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)

## Write list of filtered samples
filteredsample = filteredsample[,1:2]
excludesample = persampletable[!persampletable$SAMPLE %in% filteredsample$SAMPLE,]
excludesample = excludesample[,1:2]
write.table(filteredsample, file='filtered_sample_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(excludesample, file='exclude_sample_0805.txt', sep="\t", row.names=FALSE, quote = FALSE)

## List of samples for PheWAS
allhomofiltered = allhomo[allhomo$SAMPLE %in% filteredsample$SAMPLE,]
int_variant = c('3243_A_G','8344_A_G','8356_T_C','8363_G_A','3460_G_A','11778_G_A','14484_T_C')
int_sample_het = allhetfiltered[allhetfiltered$unique %in% int_variant,]
int_sample_homo = allhomofiltered[allhomofiltered$unique %in% int_variant,]
int_sample = rbind(int_sample_het, int_sample_homo)
int_sample = transform(int_sample, identifier=as.numeric(factor(int_sample$unique)))
int_sample = int_sample[,c(5:30,32:35)]
write.table(int_sample, file='pathogenic_variants_0809.txt', sep="\t", row.names=FALSE, quote = FALSE)

## Read in
perpostable = fread('per_pos_ref_alt_count_with_all_filter_0805.txt')
allhet = fread('all_variants_het_allfilter_0805.txt')
allhomo = fread('all_variants_homo_allfilter_0805.txt')
persampletable = fread('per_sample_count_with_all_filter_0805.txt')
