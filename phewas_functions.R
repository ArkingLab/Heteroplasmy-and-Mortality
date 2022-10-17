# Functions for PheWAS analysis
# Create analysis data frame for passing into process.glm()
# Input is a tibble with a named column `IID` = a list of variant carrier ids
create_phewas_data.input <- function(var.carriers, subset.to.unrelated=TRUE, subset.to.white.British.ancestry=TRUE, filter.to.exomes=TRUE,
                                     ukb.center.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/ukb.center.rds',
                                     ukb.icd10_phecodes.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.icd10_phecodes.rds',
                                     ukb.dat.path='/dcl01/arking/data/active/ukb_mitoscore_pheWAS/n502485.ukb.data.rds'){
  # Helper function
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
  
  # read pieces of data to form dataset
  message('reading pre-req data')
  ukb.center <- getobj(ukb.center.path)
  i <- getobj(ukb.icd10_phecodes.path) # phecodes
  dat <- getobj(ukb.dat.path)$ukb.data # ukb data
  
  # ukb.center <- getobj('~/ukb/UKB -mtDNA_2022/data_tables/ukb.center.rds')
  # i <- getobj('~/ukb/UKB -mtDNA_2022/data_tables//n502485.icd10_phecodes.rds') # phecodes
  # dat <- getobj('~/ukb/UKB -mtDNA_2022/data_tables/n502485.ukb.data.rds')$ukb.data # ukb data
  
  # ukb.center <- getobj('~/mtrv/ukb.center.rds')
  # i <- getobj('~/mtrv/n502485.icd10_phecodes.rds') # phecodes
  # dat <- getobj('~/mtrv/n502485.ukb.data.rds')$ukb.data # ukb data
  
  # create synthetic allele of rv carriers
  message('creating synthetic allele from var.carriers')
  var <- rep(0, nrow(dat))
  var[which(dat$id %in% var.carriers$IID)] <- 1
  var <- tibble(id=dat$id, var=as.factor(var))
  
  # create dataset
  message('creating dataset')
  data <- dat %>% 
    inner_join(i, by='id') %>% 
    inner_join(ukb.center, by=c('id'='IID')) %>% 
    inner_join(var, by='id') %>% 
    dplyr::select(id,age, sex, Center, var, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation, 
                  all_of(paste0('PC',1:40,separate='')), 
                  all_of(colnames(i)[2:ncol(i)]))
  data$Center <- as.factor(data$Center)
  data$sex <- as.factor(data$sex)
  data$genotyping.array <- as.factor(data$genotyping.array)
  data$in.white.British.ancestry.subset <- as.factor(data$in.white.British.ancestry.subset)
  data$used.in.pca.calculation <- as.factor(data$used.in.pca.calculation)
  data$age2 <- data$age^2
  
  # # filter datset to exome IIDs only
  # # dx download vkp/450k_exomes.IIDs.list
  # message('filtering dataset')
  # exomes.iids <- readLines('~/mtrv/450k_exomes.IIDs.list')
  # data <- data %>% filter(id %in% exomes.iids)
  
  # filter datset to exome IIDs only
  if(filter.to.exomes){
    # dx download vkp/450k_exomes.IIDs.list
    message('filtering dataset')
    exomes.iids <- readLines('~/mtrv/450k_exomes.IIDs.list')
    data <- data %>% filter(id %in% exomes.iids)
  }
  
  message('number of var: ')
  data %>% count(var) %>% print()
  
  # subset to unrelated IIDs only (optional)
  if(subset.to.unrelated){
    data <- data %>% filter(used.in.pca.calculation == 1)
    message('number of var: ')
    data %>% count(var) %>% print()
  }
  
  # subset to white IIDs only (optional)
  if(subset.to.white.British.ancestry) {
    data <- data %>% filter(in.white.British.ancestry.subset == 1)
    message('number of var: ')
    data %>% count(var)
  }
  
  # Compute sex.check
  message('computing additional data, and filtering phecodes')
  sex.check <-
    apply(data %>% dplyr::select(all_of(colnames(i)[2:ncol(i)])), 2, function(x) {
      ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
    })
  
  # Compute n.var.cases to figure out which phecodes to ignore for testing
  tmp_data <- data %>% select(var, all_of(colnames(i)[2:ncol(i)]))
  var_tmp <- tmp_data$var
  n.var.cases <- apply(tmp_data[,2:ncol(tmp_data)], 2, 
                       # function(x) length(which(x==TRUE & var==1)))
                       function(x) length(which(x==TRUE & var_tmp==1)))
  # n.var.cases <- n.var.cases[-1]
  
  # Get list of phecodes and 
  # filter phecodes such that at least one variant carrier is a case
  # i.e. number of cases with variant > 0 (reduces number of tests)
  phecodes <- colnames(i)[2:ncol(i)]
  phecodes <- phecodes[-which(n.var.cases == 0)]
  # phecodes <- phecodes[-which(n.cases <= 20)] # alt filter
  
  return(list(data, sex.check, phecodes))
}

# Create analysis data frame for passing into process.glm()
# Input is a tibble with a column `IID` 
# and a corresponding scored synthetic allele column
create_phewas_data.input2 <- function(var.carriers){
  # Helper function
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
  
  # read pieces of data to form dataset
  message('reading pre-req data')
  ukb.center <- getobj('~/ukb/UKB -mtDNA_2022/data_tables/ukb.center.rds')
  i <- getobj('~/ukb/UKB -mtDNA_2022/data_tables//n502485.icd10_phecodes.rds') # phecodes
  dat <- getobj('~/ukb/UKB -mtDNA_2022/data_tables/n502485.ukb.data.rds')$ukb.data # ukb data
  
  # ukb.center <- getobj('~/mtrv/ukb.center.rds')
  # i <- getobj('~/mtrv/n502485.icd10_phecodes.rds') # phecodes
  # dat <- getobj('~/mtrv/n502485.ukb.data.rds')$ukb.data # ukb data
  
  # create synthetic allele of rv carriers
  message('attaching scored allele from var.carriers')
  var <- rep(0, nrow(dat))
  # var[which(dat$id %in% var.carriers$IID)] <- 1
  var[which(dat$id %in% var.carriers$IID)] <- var.carriers$scored.allele
  # var <- tibble(id=dat$id, var=as.factor(var))
  var <- tibble(id=dat$id, var=var)
  
  # create dataset
  message('creating dataset')
  data <- dat %>% 
    inner_join(i, by='id') %>% 
    inner_join(ukb.center, by=c('id'='IID')) %>% 
    inner_join(var, by='id') %>% 
    dplyr::select(id,age, sex, Center, var, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation, 
                  all_of(paste0('PC',1:40,separate='')), 
                  all_of(colnames(i)[2:ncol(i)]))
  data$Center <- as.factor(data$Center)
  data$sex <- as.factor(data$sex)
  data$genotyping.array <- as.factor(data$genotyping.array)
  data$in.white.British.ancestry.subset <- as.factor(data$in.white.British.ancestry.subset)
  data$used.in.pca.calculation <- as.factor(data$used.in.pca.calculation)
  data$age2 <- data$age^2
  
  # filter datset to exome IIDs only
  # dx download vkp/450k_exomes.IIDs.list
  message('filtering dataset')
  exomes.iids <- readLines('~/mtrv/450k_exomes.IIDs.list')
  data <- data %>% filter(id %in% exomes.iids)
  
  # subset to white unrelated IIDs only
  data <- data %>% filter(in.white.British.ancestry.subset==1 & used.in.pca.calculation==1)
  
  # Compute sex.check
  message('computing additional data, and filtering phecodes')
  sex.check <-
    apply(data %>% dplyr::select(all_of(colnames(i)[2:ncol(i)])), 2, function(x) {
      ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
    })
  
  # Compute n.var.cases to figure out which phecodes to ignore for testing
  tmp_data <- data %>% select(var, all_of(colnames(i)[2:ncol(i)]))
  var_tmp <- tmp_data$var
  n.var.cases <- apply(tmp_data[,2:ncol(tmp_data)], 2,
                       ## function(x) length(which(x==TRUE & var==1)))
                       function(x) length(which(x==TRUE & var_tmp!=0))) # <--- changed from create_phewas..input()
  ## n.var.cases <- n.var.cases[-1]

  # Get list of phecodes and 
  # filter phecodes such that at least one variant carrier is a case
  # i.e. number of cases with variant > 0 (reduces number of tests)
  phecodes <- colnames(i)[2:ncol(i)]
  phecodes <- phecodes[-which(n.var.cases == 0)]
  # phecodes <- phecodes[-which(n.cases <= 20)] # alt filter
  
  return(list(data, sex.check, phecodes))
}


# Create analysis data frame for passing into process.glm()
# Input is a tibble with a column `IID` and
# a corresponding 'qvar' column that contains a quantitative phenotype of 
# interest -- i.e. this is the independent variable in a model of
# phecode ~ qvar + cov
create_phewas_data.input3 <- function(var.carriers, filter.to.exomes=FALSE, subset.to.unrelated=FALSE, subset.to.white.British.ancestry=FALSE,
                                      ukb.center.path='~/ukb/UKB -mtDNA_2022/data_tables/ukb.center.rds',
                                      ukb.icd10_phecodes.path='~/ukb/UKB -mtDNA_2022/data_tables//n502485.icd10_phecodes.rds',
                                      ukb.dat.path='~/ukb/UKB -mtDNA_2022/data_tables/n502485.ukb.data.rds'){
  # Helper function
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
  
  # read pieces of data to form dataset
  message('reading pre-req data')
  ukb.center <- getobj(ukb.center.path)
  i <- getobj(ukb.icd10_phecodes.path) # phecodes
  dat <- getobj(ukb.dat.path)$ukb.data # ukb data
  
  # ukb.center <- getobj('~/mtrv/ukb.center.rds')
  # i <- getobj('~/mtrv/n502485.icd10_phecodes.rds') # phecodes
  # dat <- getobj('~/mtrv/n502485.ukb.data.rds')$ukb.data # ukb data
  
  # create synthetic allele of rv carriers
  message('attaching qvar from var.carriers')
  # var <- rep(0, nrow(dat))
  # # var[which(dat$id %in% var.carriers$IID)] <- 1
  # var[which(dat$id %in% as.double(var.carriers$IID))] <- var.carriers$qvar
  # # var <- tibble(id=dat$id, var=as.factor(var))
  # var <- tibble(id=dat$id, var=var)

  var <- tibble(id=dat$id)
  var <- var %>% left_join(var.carriers, by=c('id'='IID')) %>% dplyr::select(IID=id, var=qvar)

    
  # create dataset
  message('creating dataset')
  data <- dat %>% 
    inner_join(i, by='id') %>% 
    inner_join(ukb.center, by=c('id'='IID')) %>% 
    inner_join(var, by=c('id'='IID')) %>% 
    dplyr::select(id,age, sex, Center, var, genotyping.array, in.white.British.ancestry.subset, used.in.pca.calculation, 
                  all_of(paste0('PC',1:40,separate='')), 
                  all_of(colnames(i)[2:ncol(i)]))
  data$Center <- as.factor(data$Center)
  data$sex <- as.factor(data$sex)
  data$genotyping.array <- as.factor(data$genotyping.array)
  data$in.white.British.ancestry.subset <- as.factor(data$in.white.British.ancestry.subset)
  data$used.in.pca.calculation <- as.factor(data$used.in.pca.calculation)
  data$age2 <- data$age^2
  
  # filter datset to exome IIDs only
  if(filter.to.exomes){
    # dx download vkp/450k_exomes.IIDs.list
    message('filtering dataset')
    exomes.iids <- readLines('~/mtrv/450k_exomes.IIDs.list')
    data <- data %>% filter(id %in% exomes.iids)
  }
  
  # subset to unrelated IIDs only (optional)
  if(subset.to.unrelated){
    data <- data %>% filter(used.in.pca.calculation == 1)
    message('number of var: ')
    data %>% count(var!=0) %>% print()
  }
  
  # subset to white IIDs only (optional)
  if(subset.to.white.British.ancestry) {
    data <- data %>% filter(in.white.British.ancestry.subset == 1)
    message('number of var: ')
    data %>% count(var!=0)
  }
  
  # Compute sex.check
  message('computing additional data, and filtering phecodes')
  sex.check <-
    apply(data %>% dplyr::select(all_of(colnames(i)[2:ncol(i)])), 2, function(x) {
      ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE) # set true if both sexes present
    })
  
  # Compute n.var.cases to figure out which phecodes to ignore for testing
  tmp_data <- data %>% select(var, all_of(colnames(i)[2:ncol(i)]))
  var_tmp <- tmp_data$var
  n.var.cases <- apply(tmp_data[,2:ncol(tmp_data)], 2,
                       ## function(x) length(which(x==TRUE & var==1)))
                       function(x) length(which(x==TRUE & var_tmp!=0))) # <--- changed from create_phewas..input()
  ## n.var.cases <- n.var.cases[-1]
  
  # Get list of phecodes and 
  # filter phecodes such that at least one variant carrier is a case
  # i.e. number of cases with variant > 0 (reduces number of tests)
  phecodes <- colnames(i)[2:ncol(i)]
  phecodes <- phecodes[-which(n.var.cases == 0)]
  # phecodes <- phecodes[-which(n.cases <= 20)] # alt filter
  
  return(list(data, sex.check, phecodes))
}

# Run PheWAS on synthetic allele
process.glm <- function(phecodes, both_sexes=NULL, data_desc, ncores=1, var='var', 
                        file.conn=NULL,back_correct_FirthSE=TRUE, debug_log=F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    header=c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls',
             'n.var.controls','case.ctrl_ratio','firth_correct')
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # ## Create a bigmatrix object from `data`
  # data2 <- bigmemory::as.big.matrix(as.data.frame(data))
  # data_desc <- bigmemory::describe(data2)
  # rm(data)
  # doParallel::registerDoParallel(ncores)
  # ifelse(debug_log,
  #        cl = parallel::makeCluster(ncores, outfile = ''),
  #        cl = parallel::makeCluster(ncores))
  cl = parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  
  x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
               .noexport = c('sex.check'),
               # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
               # .combine = 'rbind',
               .packages = c('splines','logistf','biganalytics')) %dopar% {
                 data_ <- bigmemory::attach.big.matrix(data_desc)
                 
                 phecode <- phecodes[i]
                 if(is.null(both_sexes)){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else if(!both_sexes[as.character(phecode)]){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=2) + Center + ',
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else {
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 }
                 message(f)
                 
                 # Run GLM
                 fit <-
                   biganalytics::bigglm.big.matrix(
                     formula(f),
                     data_,
                     family=binomial(link='logit'),
                     fc=ifelse(!is.null(both_sexes), 
                               yes = ifelse(!both_sexes[as.character(phecode)], 
                                      c('var','Center'), 
                                      c('var','sex','Center')),
                               no = c('var','sex','Center')),
                     quiet=T
                   )
                 # Firth Correct GLM p-value if it is <= 0.05
                 firth_correct=FALSE
                 # if(summary(fit)$coef[2,4] <= 0.05){
                 if(summary(fit)$mat[2,'p'] <= 0.05){
                   message(phecode, ' Firth Correcting...')
                   firth_correct=TRUE
                   fit.firth <-
                     logistf(
                       formula(f),
                       as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
                       control = logistf.control(maxit = 100, maxstep = 50),
                       plcontrol = logistpl.control(maxit = 400, maxstep = 20),
                       plconf = c(2),
                       model=F
                     )
                 }
                 
                 # extract values from GLM or Firth model -> write to file, and/or return
                 val <- rep(NA, 12)
                 names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases',
                                 'n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
                 val['phecode'] <- phecode
                 if(firth_correct) {
                   val['beta'] <- fit.firth$coefficients[2]
                   val['OR'] <- exp(fit.firth$coefficients[2])
                   val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
                   val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
                   val['p'] <- fit.firth$prob[2]
                   if(back_correct_FirthSE){
                     val['SE'] <- # recompute SE using beta and Chisq value
                       abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
                   }
                 } else {
                   # val['beta'] <- coef(fit)[2]
                   # val['OR'] <- exp(coef(fit)[2])
                   # val['SE'] <- summary(fit)$coef[2, 2]
                   # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
                   # val['p'] <- summary(fit)$coef[2, 4]
                   val['beta'] <- summary(fit)$mat[2,'Coef']
                   val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
                   val['SE'] <- summary(fit)$mat[2,'SE']
                   val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
                   val['p'] <- summary(fit)$mat[2,'p']
                 }
                 # val['n.cases'] <- n.cases[phecode]
                 # data <- data %>% select(var, all_of(phecode))
                 # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
                 # val['n.controls'] <- n.controls[phecode]
                 # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
                 # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
                 # val['firth_correct'] <- firth_correct
                 # message('table var: ', unique(data_[,'var']))
                 val['n.cases'] <- length(which(data_[,phecode]==1)) # n.cases[phecode]
                 # val['n.var.cases'] <- length(which(data_[,'var']==1 & data_[, phecode]==1)) # note: `var`==1 indicates presence of variant
                 val['n.var.cases'] <- length(which(data_[,'var']==2 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant
                 
                 val['n.controls'] <- length(which(data_[,phecode]==0))
                 # val['n.var.controls'] <- length(which(data_[,'var']==1 & data_[,phecode]==0)) # note: `var`==1 indicates presence of variant
                 val['n.var.controls'] <- length(which(data_[,'var']==2 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant
                 
                 val['case.ctrl_ratio'] <- 
                   (length(which(data_[,phecode]==1))) / (length(which(data_[,phecode]==0)))
                 val['firth_correct'] <- firth_correct
                 
                 # optionally, stream to file
                 if(!is.null(file.conn)){
                   val_ <- paste(c(val), sep='', collapse='\t')
                   readr::write_lines(x = val_, file = file.conn, append = TRUE)
                 }
                 return(val)
               }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  closeAllConnections()
  x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Run PheWAS on scored synthetic allele
process.glm2 <- function(phecodes, both_sexes=NULL, data_desc, ncores=1, var='var', 
                        file.conn=NULL,back_correct_FirthSE=TRUE, debug_log=F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    header=c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls',
             'n.var.controls','case.ctrl_ratio','firth_correct')
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # ## Create a bigmatrix object from `data`
  # data2 <- bigmemory::as.big.matrix(as.data.frame(data))
  # data_desc <- bigmemory::describe(data2)
  # rm(data)
  # doParallel::registerDoParallel(ncores)
  # ifelse(debug_log,
  #        cl = parallel::makeCluster(ncores, outfile = ''),
  #        cl = parallel::makeCluster(ncores))
  cl = parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  
  x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
               .noexport = c('sex.check'),
               # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
               # .combine = 'rbind',
               .packages = c('splines','logistf','biganalytics')) %dopar% {
                 data_ <- bigmemory::attach.big.matrix(data_desc)
                 
                 phecode <- phecodes[i]
                 if(is.null(both_sexes)){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else if(!both_sexes[as.character(phecode)]){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=2) + Center + ',
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else {
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 }
                 message(f)
                 
                 # Run GLM
                 fit <-
                   biganalytics::bigglm.big.matrix(
                     formula(f),
                     data_,
                     family=gaussian(link='identity'),
                     fc=ifelse(!is.null(both_sexes), 
                               yes = ifelse(!both_sexes[as.character(phecode)], # <--- DIFF
                                            c('Center'),  # <--- DIFF
                                            c('sex','Center')), # <--- DIFF
                               no = c('sex','Center')), # <--- DIFF
                     quiet=T
                   )
                 # Firth Correct GLM p-value if it is <= 0.05
                 firth_correct=FALSE
                 # if(summary(fit)$coef[2,4] <= 0.05){
                 if(summary(fit)$mat[2,'p'] <= 0.05){
                   message(phecode, ' Firth Correcting...')
                   firth_correct=TRUE
                   fit.firth <-
                     logistf(
                       formula(f),
                       as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
                       control = logistf.control(maxit = 100, maxstep = 50),
                       plcontrol = logistpl.control(maxit = 400, maxstep = 20),
                       plconf = c(2),
                       model=F
                     )
                 }
                 
                 # extract values from GLM or Firth model -> write to file, and/or return
                 val <- rep(NA, 12)
                 names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases',
                                 'n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
                 val['phecode'] <- phecode
                 if(firth_correct) {
                   val['beta'] <- fit.firth$coefficients[2]
                   val['OR'] <- exp(fit.firth$coefficients[2])
                   val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
                   val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
                   val['p'] <- fit.firth$prob[2]
                   if(back_correct_FirthSE){
                     val['SE'] <- # recompute SE using beta and Chisq value
                       abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
                   }
                 } else {
                   # val['beta'] <- coef(fit)[2]
                   # val['OR'] <- exp(coef(fit)[2])
                   # val['SE'] <- summary(fit)$coef[2, 2]
                   # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
                   # val['p'] <- summary(fit)$coef[2, 4]
                   val['beta'] <- summary(fit)$mat[2,'Coef']
                   val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
                   val['SE'] <- summary(fit)$mat[2,'SE']
                   val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
                   val['p'] <- summary(fit)$mat[2,'p']
                 }
                 # val['n.cases'] <- n.cases[phecode]
                 # data <- data %>% select(var, all_of(phecode))
                 # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
                 # val['n.controls'] <- n.controls[phecode]
                 # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
                 # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
                 # val['firth_correct'] <- firth_correct
                 # message('table var: ', unique(data_[,'var']))
                 val['n.cases'] <- length(which(data_[,phecode]==1)) # n.cases[phecode]
                 # val['n.var.cases'] <- length(which(data_[,'var']==1 & data_[, phecode]==1)) # note: `var`==1 indicates presence of variant
                 # val['n.var.cases'] <- length(which(data_[,'var']==2 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant
                 # val['n.var.cases'] <- length(which(data_[,'var']!=0 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
                 val['n.var.cases'] <- mean(data_$var[which(data_[,'var']!=0 & data_[, phecode]==1)], na.rm=T) # <--- DIFF from process.glm()
                 
                 val['n.controls'] <- length(which(data_[,phecode]==0))
                 # val['n.var.controls'] <- length(which(data_[,'var']==1 & data_[,phecode]==0)) # note: `var`==1 indicates presence of variant
                 # val['n.var.controls'] <- length(which(data_[,'var']==2 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant
                 val['n.var.controls'] <- length(which(data_[,'var']!=0 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
                 
                 val['case.ctrl_ratio'] <- 
                   (length(which(data_[,phecode]==1))) / (length(which(data_[,phecode]==0)))
                 val['firth_correct'] <- firth_correct
                 
                 # optionally, stream to file
                 if(!is.null(file.conn)){
                   val_ <- paste(c(val), sep='', collapse='\t')
                   readr::write_lines(x = val_, file = file.conn, append = TRUE)
                 }
                 return(val)
               }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  closeAllConnections()
  x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Run PheWAS on scored synthetic allele # NEWEST (as of 6/25/22)
process.glm2_rvscore <- function(phecodes, both_sexes=NULL, data_desc, ncores=1, var='var', 
                                 file.conn=NULL,back_correct_FirthSE=TRUE, 
                                 var.is_factor=FALSE, debug_log=F) {
  require(dplyr)
  require(logistf)
  # require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    header=c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls',
             'n.var.controls','case.ctrl_ratio','firth_correct')
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # ## Create a bigmatrix object from `data`
  # data2 <- bigmemory::as.big.matrix(as.data.frame(data))
  # data_desc <- bigmemory::describe(data2)
  # rm(data)
  # doParallel::registerDoParallel(ncores)
  # ifelse(debug_log,
  #        cl = parallel::makeCluster(ncores, outfile = ''),
  #        cl = parallel::makeCluster(ncores))
  cl = parallel::makeCluster(ncores, outfile="")
  doSNOW::registerDoSNOW(cl)
  # registerDoParallel(cl)
  message('ncores: ', ncores)
  
  x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
               .noexport = c('sex.check'),
               # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
               # .combine = 'rbind',
               .packages = c('splines','logistf','biganalytics')) %dopar% {
                 data_ <- bigmemory::attach.big.matrix(data_desc)
                 
                 phecode <- phecodes[i]
                 if(is.null(both_sexes)){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else if(!both_sexes[as.character(phecode)]){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=2) + Center + ',
                            paste0('PC',1:40,separate='',collapse='+'))
                 } else {
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
                            paste0('PC',1:40,separate='',collapse='+'))
                 }
                 message(f)
                 # cat(f)
                 
                 # vector of factor variable names for bigglm.big.matrix
                 var.is_factor=T
                 get.fc <- function(both_sexes, phecode, var.is_factor) {
                   fc = NULL
                   if (!both_sexes[as.character(phecode)]) {
                     fc = c('Center')
                   } else {
                     fc = c('sex', 'Center')
                   }
                   if(var.is_factor){
                     fc <- c(fc, 'var')
                   }
                   return(fc)
                 }
                 
                 # Run GLM
                 fit <-
                   biganalytics::bigglm.big.matrix(
                     formula(f),
                     data_,
                     # family=gaussian(link='identity'),
                     family=binomial(link='logit'),
                     # fc=ifelse(!is.null(both_sexes), 
                     #           yes = ifelse(!both_sexes[as.character(phecode)], # <--- DIFF
                     #                        c('Center'),  # <--- DIFF
                     #                        c('sex','Center')), # <--- DIFF
                     #           no = c('sex','Center')), # <--- DIFF
                     fc=get.fc(both_sexes, phecode, var.is_factor),
                     quiet=F
                   )
                 # Firth Correct GLM p-value if it is <= 0.05
                 firth_correct=FALSE
                 # if(summary(fit)$coef[2,4] <= 0.05){
                 if(summary(fit)$mat[2,'p'] <= 0.0005){
                   message('FIRTH Correcting...', f)
                   firth_correct=TRUE
                   fit.firth <-
                     logistf(
                       formula(f),
                       as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
                       control = logistf.control(maxit = 100, maxstep = 50),
                       plcontrol = logistpl.control(maxit = 400, maxstep = 20),
                       plconf = c(2),
                       model=F
                     )
                 }
                 
                 # extract values from GLM or Firth model -> write to file, and/or return
                 val <- rep(NA, 12)
                 names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases',
                                 'n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
                 val['phecode'] <- phecode
                 if(firth_correct) {
                   val['beta'] <- fit.firth$coefficients[2]
                   val['OR'] <- exp(fit.firth$coefficients[2])
                   val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
                   val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
                   val['p'] <- fit.firth$prob[2]
                   if(back_correct_FirthSE){
                     val['SE'] <- # recompute SE using beta and Chisq value
                       abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
                   }
                 } else {
                   # val['beta'] <- coef(fit)[2]
                   # val['OR'] <- exp(coef(fit)[2])
                   # val['SE'] <- summary(fit)$coef[2, 2]
                   # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
                   # val['p'] <- summary(fit)$coef[2, 4]
                   val['beta'] <- summary(fit)$mat[2,'Coef']
                   val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
                   val['SE'] <- summary(fit)$mat[2,'SE']
                   val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
                   val['p'] <- summary(fit)$mat[2,'p']
                 }
                 # val['n.cases'] <- n.cases[phecode]
                 # data <- data %>% select(var, all_of(phecode))
                 # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
                 # val['n.controls'] <- n.controls[phecode]
                 # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
                 # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
                 # val['firth_correct'] <- firth_correct
                 # message('table var: ', unique(data_[,'var']))
                 val['n.cases'] <- length(which(data_[,phecode]==1)) # n.cases[phecode]
                 # val['n.var.cases'] <- length(which(data_[,'var']==1 & data_[, phecode]==1)) # note: `var`==1 indicates presence of variant
                 # val['n.var.cases'] <- length(which(data_[,'var']==2 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant
                 # val['n.var.cases'] <- length(which(data_[,'var']!=0 & data_[, phecode]==1)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
                 # val['n.var.cases'] <- mean(data_[,'var'][which(data_[,'var']!=0 & data_[, phecode]==1)], na.rm=T) # <--- DIFF from process.glm()
                 
                 val['n.controls'] <- length(which(data_[,phecode]==0))
                 # val['n.var.controls'] <- length(which(data_[,'var']==1 & data_[,phecode]==0)) # note: `var`==1 indicates presence of variant
                 # val['n.var.controls'] <- length(which(data_[,'var']==2 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant
                 # val['n.var.controls'] <- length(which(data_[,'var']!=0 & data_[,phecode]==0)) # note: `var`==2 indicates presence of variant <--- DIFF from process.glm()
                 
                 val['case.ctrl_ratio'] <- 
                   (length(which(data_[,phecode]==1))) / (length(which(data_[,phecode]==0)))
                 val['firth_correct'] <- firth_correct
                 
                 # optionally, stream to file
                 if(!is.null(file.conn)){
                   val_ <- paste(c(val), sep='', collapse='\t')
                   readr::write_lines(x = val_, file = file.conn, append = TRUE)
                 }
                 return(val)
               }
  
  # doParallel::stopImplicitCluster()
  stopCluster(cl)
  # closeAllConnections()
  x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Run PheWAS on a quantitative variable (used for mtscore_sum_resid) <<-- USED IN HET PAPER PHEWAS (Main + PathogenicVariant)
process.glm3 <- function(phecodes, both_sexes=NULL, data_, ncores=1, var='var', 
                         file.conn=NULL, run.firth=TRUE, back_correct_FirthSE=TRUE, 
                         var.is_factor=FALSE, missing_var=NULL, debug_log=F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    if(var.is_factor){
      header=c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases','n.controls',
               'n.var.controls','case.ctrl_ratio','firth_correct')      
    } else {
      header=c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases','n.controls',
               'mean.var.controls','case.ctrl_ratio','firth_correct')
    }
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  # data_ <- bigmemory::attach.big.matrix(data_desc)
  for(i in 1:length(phecodes)){
                 phecode <- phecodes[i]
                 
                 # # Formula with het_count
                 # if(is.null(both_sexes)){
                 #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + het_count')#,
                 #            # paste0('PC',1:40,separate='',collapse='+'))
                 # } else if(!both_sexes[as.character(phecode)]){
                 #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center + het_count')#,
                 #            # paste0('PC',1:40,separate='',collapse='+'))
                 # } else {
                 #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + het_count')#,
                 #            # paste0('PC',1:40,separate='',collapse='+'))
                 # }
                 
                 # Formula without het_count
                 if(is.null(both_sexes)){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + SmokingStatus')#,
                   # paste0('PC',1:40,separate='',collapse='+'))
                 } else if(!both_sexes[as.character(phecode)]){
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center + SmokingStatus')#,
                   # paste0('PC',1:40,separate='',collapse='+'))
                 } else {
                   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center + SmokingStatus')#,
                   # paste0('PC',1:40,separate='',collapse='+'))
                 }
                 message(f)
                 
                 # Run GLM
                 dat=as.data.frame(data_[, c('id', phecode, 'var','age','sex','Center','het_count','SmokingStatus')])
                 dat$sex <- as.factor(dat$sex)
                 dat$Center <- as.factor(dat$Center)
                 dat$het_count <- as.factor(dat$het_count)
                 dat$SmokingStatus <- as.factor(dat$SmokingStatus)
                 if(!is.null(missing_var)){ # set missing var to NA if IIDs given
                   dat$var[which(dat$id %in% missing_var)] <- NA
                 }
                 if(var.is_factor){
                   dat$var <- as.factor(dat$var)
                 }
                
                 fit <- 
                   glm(formula = formula(f), 
                       data=dat,
                       family=binomial(link='logit'))
                 
                 # print(summary(fit))

                 # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
                 # fc <- c('sex','Center')
                 # if(!is.null(both_sexes)) {
                 #   if (!both_sexes[as.character(phecode)]) {
                 #     fc <- fc[!fc %in% c('sex')]
                 #   }
                 # }
                 # fit <-
                 #   biganalytics::bigglm.big.matrix(
                 #     formula(f),
                 #     data_,
                 #     family=gaussian(link='identity'),
                 #     fc=fc,
                 #     # fc=ifelse(!is.null(both_sexes), 
                 #     #           yes = ifelse(!both_sexes[as.character(phecode)],
                 #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
                 #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
                 #     #           no = c('sex','Center')),
                 #     quiet=T
                 #   )
                 
                 # Firth Correct GLM p-value if it is <= 0.05
                 firth_correct=FALSE
                 if(run.firth){
                 # if(summary(fit)$coef[2,4] <= 0.05){
                 if(summary(fit)$coef[2,4] <= 5e-02){
                 # if(summary(fit)$mat[2,'p'] <= 0.05){
                   message(phecode, ' Firth Correcting...')
                   firth_correct=TRUE
                   fit.firth <-
                     logistf(
                       formula(f),
                       # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
                       # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
                       dat,
                       control = logistf.control(maxit = 100, maxstep = 50),
                       plcontrol = logistpl.control(maxit = 400, maxstep = 20),
                       plconf = c(2),
                       model=F
                     )
                 }
                 }
                 
                 # extract values from GLM or Firth model -> write to file, and/or return
                 val <- rep(NA, 12)
                 if(var.is_factor){
                   names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','n.var.cases',
                                   'n.controls','n.var.controls','case.ctrl_ratio','firth_correct')
                 } else {
                   names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases',
                                   'n.controls','mean.var.controls','case.ctrl_ratio','firth_correct')
                 }
                 val['phecode'] <- phecode
                 if(firth_correct) {
                   val['beta'] <- fit.firth$coefficients[2]
                   val['OR'] <- exp(fit.firth$coefficients[2])
                   val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
                   val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
                   val['p'] <- fit.firth$prob[2]
                   if(back_correct_FirthSE){
                     val['SE'] <- # recompute SE using beta and Chisq value
                       abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
                   }
                 } else {
                   # val['beta'] <- coef(fit)[2]
                   # val['OR'] <- exp(coef(fit)[2])
                   # val['SE'] <- summary(fit)$coef[2, 2]
                   # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
                   # val['p'] <- summary(fit)$coef[2, 4]
                   # val['beta'] <- summary(fit)$mat[2,'Coef']
                   val['beta'] <- summary(fit)$coef[2,'Estimate']
                   # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
                   val['OR'] <- exp(summary(fit)$coef[2,'Estimate'])
                   # val['SE'] <- summary(fit)$mat[2,'SE']
                   val['SE'] <- summary(fit)$coef[2,'Std. Error']
                   # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
                   val['chisq'] <- qchisq(1 - summary(fit)$coef[2,'Pr(>|z|)'], df=1)
                   # val['p'] <- summary(fit)$mat[2,'p']
                   val['p'] <- summary(fit)$coef[2,'Pr(>|z|)']
                 }
                 # val['n.cases'] <- n.cases[phecode]
                 # data <- data %>% select(var, all_of(phecode))
                 # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
                 # val['n.controls'] <- n.controls[phecode]
                 # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
                 # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
                 # val['firth_correct'] <- firth_correct
                 # message('table var: ', unique(data_[,'var']))
                 val['n.cases'] <- length(which(data_[,phecode]==TRUE & !is.na(data_[,'var'])))
                 val['n.controls'] <- length(which(data_[,phecode]==FALSE & !is.na(data_[,'var'])))
                 if(var.is_factor){
                   val['n.var.cases'] <- length(which(data_[,phecode]==TRUE & data_[,'var']==1))
                   val['n.var.controls'] <- length(which(data_[,phecode]==FALSE & data_[,'var']==1))
                 } else {
                   var_ <- dat[,'var']
                   val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==TRUE)], na.rm=T)
                   val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==FALSE)], na.rm=T)
                 }
                 
                 val['case.ctrl_ratio'] <- 
                   (length(which(data_[,phecode]==TRUE & !is.na(data_[,'var'])))) / (length(which(data_[,phecode]==FALSE & !is.na(data_[,'var']))))
                 val['firth_correct'] <- firth_correct
                 
                 # optionally, stream to file
                 if(!is.null(file.conn)){
                   val_ <- paste(c(val), sep='', collapse='\t')
                   readr::write_lines(x = val_, file = file.conn, append = TRUE)
                 }
                 rm(dat)
                 # return(val)
                 x[[i]] <- val
               }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Run PheWAS on a quantitative variable
process.glm3.1 <- function(phecodes, both_sexes=NULL, data_desc, ncores=1, var='var', 
                         file.conn=NULL,back_correct_FirthSE=TRUE, missing_var=NULL, debug_log=F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    header=c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases','n.controls',
             'mean.var.controls','case.ctrl_ratio','firth_correct')
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  data_ <- bigmemory::attach.big.matrix(data_desc)
  for(i in 1:length(phecodes)){
    phecode <- phecodes[i]
    if(is.null(both_sexes)){
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#, 
      # paste0('PC',1:40,separate='',collapse='+'))
    } else if(!both_sexes[as.character(phecode)]){
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center')#,
      # paste0('PC',1:40,separate='',collapse='+'))
    } else {
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#, 
      # paste0('PC',1:40,separate='',collapse='+'))
    }
    message(f)
    
    # Run GLM
    dat=as.data.frame(data_[, c('id', phecode, 'var','age','sex','Center','het_count')])
    dat$sex <- as.factor(dat$sex)
    dat$Center <- as.factor(dat$Center)
    # dat$het_count <- as.factor(dat$het_count)
    if(!is.null(missing_var)){ # set missing var to NA if IIDs given
      dat$var[which(dat$id %in% missing_var)] <- NA
    }
    
    fit <- 
      glm(formula = formula(f), 
          data=dat,
          family=binomial(link='logit'))
    
    # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
    # fc <- c('sex','Center')
    # if(!is.null(both_sexes)) {
    #   if (!both_sexes[as.character(phecode)]) {
    #     fc <- fc[!fc %in% c('sex')]
    #   }
    # }
    # fit <-
    #   biganalytics::bigglm.big.matrix(
    #     formula(f),
    #     data_,
    #     family=gaussian(link='identity'),
    #     fc=fc,
    #     # fc=ifelse(!is.null(both_sexes), 
    #     #           yes = ifelse(!both_sexes[as.character(phecode)],
    #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
    #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
    #     #           no = c('sex','Center')),
    #     quiet=T
    #   )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct=FALSE
    if(FALSE){
      # if(summary(fit)$coef[2,4] <= 0.05){
      if(summary(fit)$coef[2,4] <= 0.05){
        # if(summary(fit)$mat[2,'p'] <= 0.05){
        message(phecode, ' Firth Correcting...')
        firth_correct=TRUE
        fit.firth <-
          logistf(
            formula(f),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
            dat,
            control = logistf.control(maxit = 100, maxstep = 50),
            plcontrol = logistpl.control(maxit = 400, maxstep = 20),
            plconf = c(2),
            model=F
          )
      }}
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases',
                    'n.controls','mean.var.controls','case.ctrl_ratio','firth_correct')
    val['phecode'] <- phecode
    if(firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
      val['p'] <- fit.firth$prob[2]
      if(back_correct_FirthSE){
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      # val['beta'] <- summary(fit)$mat[2,'Coef']
      val['beta'] <- summary(fit)$coef[2,'Estimate']
      # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
      val['OR'] <- exp(summary(fit)$coef[2,'Estimate'])
      # val['SE'] <- summary(fit)$mat[2,'SE']
      val['SE'] <- summary(fit)$coef[2,'Std. Error']
      # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
      val['chisq'] <- qchisq(1 - summary(fit)$coef[2,'Pr(>|z|)'], df=1)
      # val['p'] <- summary(fit)$mat[2,'p']
      val['p'] <- summary(fit)$coef[2,'Pr(>|z|)']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    var_ <- data[,'var']
    val['n.cases'] <- length(which(data_[,phecode]==1 & !is.na(data_[,'var'])))
    val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==1)], na.rm=T)
    val['n.controls'] <- length(which(data_[,phecode]==0 & !is.na(data_[,'var'])))
    val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==0)], na.rm=T)
    
    val['case.ctrl_ratio'] <- 
      (length(which(data_[,phecode]==1 & !is.na(data_[,'var'])))) / (length(which(data_[,phecode]==0 & !is.na(data_[,'var']))))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if(!is.null(file.conn)){
      val_ <- paste(c(val), sep='', collapse='\t')
      readr::write_lines(x = val_, file = file.conn, append = TRUE)
    }
    rm(dat)
    # return(val)
    x[[i]] <- val
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}

# Run PheWAS on a quantitative variable
process.glm3.2 <- function(phecodes, both_sexes=NULL, data_, ncores=1, var='var', 
                           file.conn=NULL,back_correct_FirthSE=TRUE, missing_var=NULL, debug_log=F) {
  require(dplyr)
  require(logistf)
  require(doParallel)
  require(parallel)
  require(splines)
  require(logistf)
  require(bigmemory)
  require(biganalytics)
  require(snow)
  require(doSNOW)
  
  # pre-compute some info about dataset
  # message('Precomputing data: n.var in cases, controls, etc.')
  # sex.check <-
  #   apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) {
  #     ifelse(length(unique(data$sex[which(x)])) == 2, TRUE, FALSE)
  #   })
  # var <- data %>% select(var)
  # n.cases <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x)))
  # n.var.cases <- apply(data %>% dplyr::select(var,all_of(phecodes)), 2, function(x) length(which(x==TRUE & var==1)))
  # n.var.cases <- n.var.cases[-1]
  # n.controls <- apply(data %>% dplyr::select(all_of(phecodes)), 2, function(x) length(which(x==FALSE)))
  # case.ctrl_ratio <- n.cases / n.controls
  
  # print(str(data[,1:10]))
  # print(str(n.var.cases))
  # print(str(n.var.controls))
  
  # filter phecodes to test
  # filter: number of cases  with variant (`var`) > 0
  # phecodes <- phecodes[-which(n.var.cases == 0)]
  # message('Number of phecodes: ', length(phecodes))
  
  # if streaming results to file, write a file with a header
  if(!is.null(file.conn)){
    header=c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases','n.controls',
             'mean.var.controls','case.ctrl_ratio','firth_correct')
    readr::write_lines(x = paste(header, sep='\t', collapse='\t'), 
                       file = file.conn, append = TRUE)
  }
  
  # cl = parallel::makeCluster(ncores)
  # doSNOW::registerDoSNOW(cl)
  
  # x <- foreach(i = 1:length(phecodes), .errorhandling = 'pass', 
  #              .noexport = c('sex.check'),
  #              # .export = c('both_sexes','phecodes','data_desc','back_correct_FirthSE','file.conn','var'),
  #              # .combine = 'rbind',
  #              .packages = c('splines','logistf','biganalytics')) %dopar% {
  x <- list()
  # data_ <- bigmemory::attach.big.matrix(data_desc)
  for(i in 1:length(phecodes)){
    phecode <- phecodes[i]
    # if(is.null(both_sexes)){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#, 
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # } else if(!both_sexes[as.character(phecode)]){
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=4) + Center')#,
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # } else {
    #   f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=4) + sex + Center')#, 
    #   # paste0('PC',1:40,separate='',collapse='+'))
    # }
    if(is.null(both_sexes)){
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
               paste0('PC',1:40,separate='',collapse='+'))
    } else if(!both_sexes[as.character(phecode)]){
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age,df=2) + Center + ',
               paste0('PC',1:40,separate='',collapse='+'))
    } else {
      f=paste0('`',phecode,'`',' ~ ', var, ' + ns(age, df=2) + sex + Center + ', 
               paste0('PC',1:40,separate='',collapse='+'))
    }
    message(f)
    
    # Run GLM
    # dat=as.data.frame(data_[, c('id', phecode, 'var','age','sex','Center','het_count')])
    dat=as.data.frame(data_[, c('id', phecode, 'var','age','sex','Center',paste0('PC',1:40))])
    dat$sex <- as.factor(dat$sex)
    dat$Center <- as.factor(dat$Center)
    # dat$het_count <- as.factor(dat$het_count)
    if(!is.null(missing_var)){ # set missing var to NA if IIDs given
      dat$var[which(dat$id %in% missing_var)] <- NA
    }
    
    fit <- 
      glm(formula = formula(f), 
          data=dat,
          family=binomial(link='logit'))
    
    # alternate GLM using biganalytics (doesn't properly remove missing var values, so not using)
    # fc <- c('sex','Center')
    # if(!is.null(both_sexes)) {
    #   if (!both_sexes[as.character(phecode)]) {
    #     fc <- fc[!fc %in% c('sex')]
    #   }
    # }
    # fit <-
    #   biganalytics::bigglm.big.matrix(
    #     formula(f),
    #     data_,
    #     family=gaussian(link='identity'),
    #     fc=fc,
    #     # fc=ifelse(!is.null(both_sexes), 
    #     #           yes = ifelse(!both_sexes[as.character(phecode)],
    #     #                        yes=c('Center'),  # <--- both_sexes==TRUE, so !both_sexes==FALSE; thus only one sex present, so remove sex covar
    #     #                        no=c('sex','Center')), # <--- both_sexes==FALSE, so !both_sexes==TRUE, thus both present, so keep sex covar
    #     #           no = c('sex','Center')),
    #     quiet=T
    #   )
    
    # Firth Correct GLM p-value if it is <= 0.05
    firth_correct=FALSE
    if(FALSE){
      # if(summary(fit)$coef[2,4] <= 0.05){
      if(summary(fit)$coef[2,4] <= 0.05){
        # if(summary(fit)$mat[2,'p'] <= 0.05){
        message(phecode, ' Firth Correcting...')
        firth_correct=TRUE
        fit.firth <-
          logistf(
            formula(f),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center',paste0('PC',1:40, separate=''))]),
            # as.data.frame(data_[, c(phecode, 'var','age','sex','Center')]),
            dat,
            control = logistf.control(maxit = 100, maxstep = 50),
            plcontrol = logistpl.control(maxit = 400, maxstep = 20),
            plconf = c(2),
            model=F
          )
      }}
    
    # extract values from GLM or Firth model -> write to file, and/or return
    val <- rep(NA, 12)
    names(val) <- c('phecode','beta','SE','chisq','OR','p','n.cases','mean.var.cases',
                    'n.controls','mean.var.controls','case.ctrl_ratio','firth_correct')
    val['phecode'] <- phecode
    if(firth_correct) {
      val['beta'] <- fit.firth$coefficients[2]
      val['OR'] <- exp(fit.firth$coefficients[2])
      val['SE'] <- sqrt(diag(vcov(fit.firth))[2])
      val['chisq'] <- qchisq(1-fit.firth$prob,df=1)[2]
      val['p'] <- fit.firth$prob[2]
      if(back_correct_FirthSE){
        val['SE'] <- # recompute SE using beta and Chisq value
          abs(fit.firth$coefficients[2]) / sqrt(qchisq(1 - fit.firth$prob, df = 1)[2])
      }
    } else {
      # val['beta'] <- coef(fit)[2]
      # val['OR'] <- exp(coef(fit)[2])
      # val['SE'] <- summary(fit)$coef[2, 2]
      # val['chisq'] <- qchisq(1-summary(fit)$coef[2,4], df=1)
      # val['p'] <- summary(fit)$coef[2, 4]
      # val['beta'] <- summary(fit)$mat[2,'Coef']
      val['beta'] <- summary(fit)$coef[2,'Estimate']
      # val['OR'] <- exp(summary(fit)$mat[2,'Coef'])
      val['OR'] <- exp(summary(fit)$coef[2,'Estimate'])
      # val['SE'] <- summary(fit)$mat[2,'SE']
      val['SE'] <- summary(fit)$coef[2,'Std. Error']
      # val['chisq'] <- qchisq(1 - summary(fit)$mat[2,'p'], df=1)
      val['chisq'] <- qchisq(1 - summary(fit)$coef[2,'Pr(>|z|)'], df=1)
      # val['p'] <- summary(fit)$mat[2,'p']
      val['p'] <- summary(fit)$coef[2,'Pr(>|z|)']
    }
    # val['n.cases'] <- n.cases[phecode]
    # data <- data %>% select(var, all_of(phecode))
    # val['n.var.cases'] <- length(which(data[,1]==1 & data[,2]==TRUE))
    # val['n.controls'] <- n.controls[phecode]
    # val['n.var.controls'] <- length(which(data[,1]==1 & data[,2]==FALSE))
    # val['case.ctrl_ratio'] <- n.cases[phecode] / n.controls[phecode]
    # val['firth_correct'] <- firth_correct
    # message('table var: ', unique(data_[,'var']))
    var_ <- data[,'var']
    val['n.cases'] <- length(which(data_[,phecode]==1 & !is.na(data_[,'var'])))
    val['mean.var.cases'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==1)], na.rm=T)
    val['n.controls'] <- length(which(data_[,phecode]==0 & !is.na(data_[,'var'])))
    val['mean.var.controls'] <- mean(var_$var[which(!is.na(data_[,'var']) & data_[, phecode]==0)], na.rm=T)
    
    val['case.ctrl_ratio'] <- 
      (length(which(data_[,phecode]==1 & !is.na(data_[,'var'])))) / (length(which(data_[,phecode]==0 & !is.na(data_[,'var']))))
    val['firth_correct'] <- firth_correct
    
    # optionally, stream to file
    if(!is.null(file.conn)){
      val_ <- paste(c(val), sep='', collapse='\t')
      readr::write_lines(x = val_, file = file.conn, append = TRUE)
    }
    rm(dat)
    # return(val)
    x[[i]] <- val
  }
  
  # doParallel::stopImplicitCluster()
  # stopCluster(cl)
  # closeAllConnections()
  # x <- bind_rows(x) %>% as_tibble()
  return(x)
}


# Helper Functions (Misc)
# getobj.sizes <- function(){
#   return(sort(sapply(ls(), function(x) object.size(get(x)))))
#   sort(sapply(ls(), function(x) object.size(get(x))))
# }
# print(object.size(exomes.iids), units='auto', standard='SI') # single object's size
# list of object sizes ()
# for(i in 1:length(ls())){message(ls()[i]); print(object.size(get(ls()[i])), units='auto', standard='SI')}
