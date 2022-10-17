# code to run mitoscore heteroplasmy pheWAS

library(tools)
library(dplyr)
library(readr)
# library(parallel)
library(splines)
# library(doParallel)
library(logistf)
# library(bigmemory)
# library(biganalytics)
library(doSNOW)

# source pheWAS functions
source('/dcl01/arking/data/active/ukb_mitoscore_pheWAS/phewas_functions.R')

# read in mitoscore
x <- as_tibble(read.csv('/dcl01/arking/data/active/ukb_mitoscore_pheWAS/mtscore_sum_resid_061022.txt'))
# x <- as_tibble(read.csv('~/ukb/UKB -mtDNA_2022/data_tables/mtscore_sum_resid_061022.txt'))
var.carriers <- x %>% dplyr::select(IID=eid, qvar=mtscore)

# create dataset
data <- create_phewas_data.input3(var.carriers=var.carriers,
                                  ukb.center.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/ukb.center.rds',
                                  ukb.icd10_phecodes.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.icd10_phecodes.rds',
                                  ukb.dat.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.ukb.data.rds',
                                  subset.to.white.British.ancestry=FALSE,
                                  subset.to.unrelated=FALSE,
                                  filter.to.exomes=FALSE)
# data <- create_phewas_data.input3(var.carriers=var.carriers, 
#                                   ukb.center.path='/home/dnanexus/mtrv/ukb.center.rds', 
#                                   ukb.icd10_phecodes.path='/home/dnanexus/mtrv/n502485.icd10_phecodes.rds',
#                                   ukb.dat.path='/home/dnanexus/mtrv/n502485.ukb.data.rds',
#                                   subset.to.white.British.ancestry=FALSE,
#                                   subset.to.unrelated=FALSE,
#                                   filter.to.exomes=FALSE)
# data <- create_phewas_data.input3(var.carriers=var.carriers, 
#                                   ukb.center.path='~/ukb.center.rds', 
#                                   ukb.icd10_phecodes.path='~/n502485.icd10_phecodes.rds',
#                                   ukb.dat.path='~/n502485.ukb.data.rds',
#                                   subset.to.white.British.ancestry=FALSE,
#                                   subset.to.unrelated=FALSE,
#                                   filter.to.exomes=FALSE)
# data <- create_phewas_data.input3(var.carriers=var.carriers) # using default values for ukb.* paths

phecodes = data[[3]]
sex.check=data[[2]]
data <- data[[1]]
length(phecodes)

# Add het_count
het <- as_tibble(read.csv('/dcl01/arking/data/active/ukb_mitoscore_pheWAS/het_count_061022.txt'))
# het <- as_tibble(read.csv('~/ukb/UKB -mtDNA_2022/data_tables/het_count_061022.txt'))
data <- data %>% left_join(het, by=c('id'='eid'))

# IIDs with missing mitoscore
missing_var <- data %>% anti_join(var.carriers, by=c('id'='IID')) %>% dplyr::select(IID='id')
missing_var <- missing_var$IID

## Create a bigmatrix object from `data`
# require(bigmemory)
# require(biganalytics)
# k2 <- bigmemory::as.big.matrix(as.data.frame(data))
# data_desc <- bigmemory::describe(k2)

# Run PheWAS
# phewas <- process.glm3(phecodes=sample(phecodes), both_sexes=sex.check, 
#                        data=data_desc, ncores=1,var='var',
#                        file.conn='./phewas.mtscore2.txt')
phewas <- process.glm3(phecodes=sample(phecodes), both_sexes=sex.check, # Firth correct sig
                       data_=data, ncores=1,var='var', run.firth=TRUE,
                       file.conn='./phewas.mtscore_with_het_count.txt')

# Firth correct fails to converge for these phecodes, rerun without firth correction to get glm() outputs
phecodes_ <- c('200','204','204.12','204.1')
phewas <- process.glm3(phecodes=sample(phecodes_), both_sexes=sex.check, data_=data, ncores=1,var='var', run.firth=FALSE)

# harvest results
# old - comparing betas
# t0 <- read_tsv('~/projects/mito_het_pheWAS/phewas.mtscore.txt', skip = 2)
# t0 <- t0 %>% addPhecodeInfo() %>% arrange(p)
# t0 %>% inner_join(t1, by='phecode') %>%
#   ggplot(aes(x=beta.x, y=beta.y)) +
#   geom_point() + hrbrthemes::theme_ipsum_rc() + 
#   xlab('Beta Without Adjustment for Heteroplasmy Count') + 
#   ylab('Beta After Adjustment for Heteroplasmy Count') + 
#   labs(title='Betas before/after heteroplasmy adjustment')
# 
# t1 <- read_tsv('~/projects/mito_het_pheWAS/phewas.mtscore.sig_with_hetcount.txt')


# harvest results new 
require(PheWAS)

# Plot of PheWAS with het_count adjustment:
t0 <- read_tsv('/Users/vkp/projects/mito_het_pheWAS/phewas.mtscore_sum_resid_061022_with_het_count.txt', skip=3)
t0 <- t0 %>% filter(n.cases >= 10) %>% arrange(p)
t0 <- t0 %>% addPhecodeInfo() %>% arrange(p)
t0_ <- t0 %>% select(phenotype=phecode, p, OR)

t0.plot <- phewasManhattan(t0_, OR.size=T, OR.direction=T, annotate.size = 4.5)
t0.plot + labs(title='PheWAS of Mitoscore', 
               subtitle = 'Model: <phecode> ~ mtscore_sum_resid_061022 + age + sex + Center + Heteroplasmy Count', 
               caption = 'Age modeled as ns(age, df=4)') + theme(legend.text = element_blank())

# Plot of PheWAS without het_count:
require(PheWAS)
require(readr)
require(dplyr)
t0 <- read_tsv('/Users/vkp/projects/mito_het_pheWAS/phewas.mtscore_sum_resid_061022_no_het_count.txt', skip=3)
t0 <- t0 %>% filter(n.cases >= 10) %>% arrange(p)
t0 <- t0 %>% addPhecodeInfo() %>% arrange(p)
t0_ <- t0 %>% dplyr::select(phenotype=phecode, p, OR)
t0_ <-addPhecodeInfo(t0_, groupnums = T, groupcolors = T)

# for phewas.mtscore_sum_resid_061022_no_het_count.txt specifically,
# the p-values for Leukemia-related codes after Firth correction
# are 0, which is because `logistf()` couldn't converge on a p-val.
# So, we choose to set to `glm()` output p-values as p-values for these phecodes.
# TODO: see if Firth can be run longer to converge on correct p-value
glm_pvals <- c(1.26193934580806e-56, 3.09852721736083e-44, 3.18816790927094e-40, 5.01375848867052e-37)
t0_$p[t0_$phenotype %in% c('200','204','204.12','204.1')] <- glm_pvals

# Also, split plot up due to outlier y-values
# Manually paste both plots in powerpoint/keynote
t0.plot <- phewasManhattan(t0_,
# t0.plot <- phenotypeManhattan(t0_,
                           # suggestive.line = NA,
                           OR.size=T, 
                           OR.direction=T, 
                           # point.size = 3,
                           annotate.size = 5,
                           size.x.labels = 20,
                           sort.by.category.value = T)

# Bottom plot
bottom <- t0.plot +
  labs(title='PheWAS of Mitoscore', subtitle = 'Model: <phecode> ~ mtscore_sum_resid_061022 + age + sex + Center', caption = 'Age modeled as ns(age, df=4)') +
  theme(#legend.text = element_blank(), 
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 25)) +
  ylim(c(0,16)) +
  ylab(bquote(-log[10](italic("p-value"))))
ggsave(plot = bottom, device='tiff', units='in', height=12, width=20, path = '~/projects/mito_het_pheWAS/', filename = 'Plot_pheWAS_mtscore_sum_resid_without_het_BOTTOM.tiff')

# Top plot
top <- t0.plot + 
  labs(title='PheWAS of Mitoscore', subtitle = 'Model: <phecode> ~ mtscore_sum_resid_061022 + age + sex + Center', caption = 'Age modeled as ns(age, df=4)') +
  theme(legend.text = element_blank(), 
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 25)) +
  ylim(c(35,56)) +
  ylab(bquote(-log[10](italic("p-value"))))
ggsave(plot = top, device='tiff', units='in', height=7, width=20, path = '~/projects/mito_het_pheWAS/', filename = 'Plot_pheWAS_mtscore_sum_resid_without_het_TOP.tiff')


# Correlation plot of phecodes for top phecodes
t <- readRDS('~/projects/mito_het_pheWAS/data_phecodes_tophits.rds') # p <= 0.05/1502 (pheWAS without het count)
cor <- cor(t[,10:38], use = 'pairwise.complete.obs')
corrplot(cor, order='hclust', title = 'Correlation of significant phecodes')


# Correlation plot for individuals with leukemia only does with other top hits
getobj <- function (Rdata)
{
  require(tools)
  if (tolower(file_ext(Rdata)) == "rds") {
    return(readRDS(Rdata))
  }
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata, 
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

mtscore <- read_csv('~/projects/mito_het_pheWAS/mtscore_sum_resid_061022.txt')
phecodes <- getobj('~/projects/mito_rare-variant/resources/n502485.icd10_phecodes.rds')
leukemia_codes <- c('200','200.1','202','202.2','204','204.1','204.12','204.2','204.21','204.3')
sepsis_codes <- c('038','994','994.2')
mtscore_ <- mtscore %>% inner_join(phecodes, by=c('eid'='id')) %>%
  select(eid, all_of(leukemia_codes), all_of(sepsis_codes))

mtscore_leukemia <- mtscore_ %>% filter(`200`==TRUE | `200.1`==TRUE | `202`==TRUE | `202.2`==TRUE 
                    | `204`==TRUE | `204.1`==TRUE | `204.12`==TRUE | `204.2`==TRUE | `204.21`==TRUE | `204.3`==TRUE)

cor <- cor(mtscore_leukemia[, 2:ncol(mtscore_leukemia)], use = 'pairwise.complete.obs')

# Do Luekmia cases also have sepsis? Meaning is the sepsis signal coming from hematological cancer cases?
data_ <- data %>% filter(!is.na(var))
data_ %>% count(`200`,`038`)
# # A tibble: 9 × 3
# `200` `038`      n
# <lgl> <lgl>  <int>
#   1 FALSE FALSE 179580
# 2 FALSE TRUE    4111
# 3 FALSE NA      6854
# 4 TRUE  FALSE    530
# 5 TRUE  TRUE     170
# 6 TRUE  NA        64
# 7 NA    FALSE   1643
# 8 NA    TRUE     698
# 9 NA    NA       216


# 2x2 of Sepsis (038, rows) x Leukemia (200, cols)
#         Cases 	Controls
# Cases	    170   4111
# Controls  530   179580
fisher.test(matrix(data=c(170,530,4111, 179580), 2, 2))
# > fisher.test(matrix(data=c(170,530,4111, 179580), 2, 2))
# 
# Fisher's Exact Test for Count Data
# 
# data:  matrix(data = c(170, 530, 4111, 179580), 2, 2)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  11.68299 16.72567
# sample estimates:
# odds ratio
#   14.00602
# hema.code dummy variable
# sepsis.code dummy variable
data_ <- 
  data_ %>% 
  mutate(hema.code=ifelse(`200`==TRUE | `200.1`==TRUE | `202`==TRUE | 
                            `202.2`==TRUE | `204`==TRUE | `204.1`==TRUE | 
                            `204.12`==TRUE | `204.2`==TRUE | `204.21`==TRUE | 
                            `204.3`==TRUE, 1, 0)) %>% 
  mutate(sepsis.code=ifelse(`038`==TRUE | `994`==TRUE | `994.2`==TRUE, 1, 0))

data_ %>% count(hema.code, sepsis.code)
# # A tibble: 9 × 3
# hema.code sepsis.code      n
# <dbl>       <dbl>  <int>
#   1         0           0 177545
# 2         0           1   3800
# 3         0          NA   6606
# 4         1           0   2136
# 5         1           1    862
# 6         1          NA    279
# 7        NA           0   2072
# 8        NA           1    317
# 9        NA          NA    249
fisher.test(matrix(data=c(862,2136,3800,177545), 2, 2))
# Fisher's Exact Test for Count Data
# data:  matrix(data = c(862, 2136, 3800, 177545), 2, 2)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  17.28351 20.56079
# sample estimates:
# odds ratio
#   18.85441

# Dan: "What about if you exclude hematological cancers...do we still see an association with sepsis?"
leukemia_codes <- c('200','200.1','202','202.2','204','204.1','204.12','204.2','204.21','204.3')
sepsis_codes <- c('038','994','994.2')

f <- formula('sepsis.code ~ var + ns(age, df=4) + sex + Center')
summary(glm(f, data = data_ %>% filter(hema.code==0), family=binomial(link='logit')))
# results: https://arkinglab.slack.com/archives/D50R63VGV/p1660077657400879
# results indicate the sepsis signal is still significant (p=0.00392) even after removing hema.codes
# further adjusting for het_count, this p=0.0539)

# PHESANT analysis (run by Dan)
g <- read_tsv('/Users/vkp/Library/Mobile Documents/com~apple~CloudDocs/Downloads/phesant_results_combined_072622.txt')
g %>% mutate(bonf.sig=ifelse(pvalue <= (0.05/nrow(g)), 1, 0)) %>% filter(bonf.sig==1) %>% View()


# PHESANT analysis visualization


