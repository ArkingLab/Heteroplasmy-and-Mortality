ibrary(RColorBrewer)
library(pheatmap)
complex<-c("ALL","DLOOP","RRNA","TRNA","ComplexI","ComplexIII","ComplexIV","ComplexV")
folder<-"C:/Users/darking/OneDrive - Johns Hopkins/darking/papers/current/mito.het_UKB/results/"
##more memory efficient
##remove large negative value betas (wonky stats due to small samples)
filename<-paste(folder,"phesant_",complex[1],".txt",sep="")
phesant<-read_tsv(file=filename)%>% select(-(isTraitOfInterest:Path),-(lower:upper)) %>% filter(beta> -2)
names(phesant)[4:5]<-paste0(complex[1],c("_beta","_pvalue"))
for(i in 2:length(complex)) {
  name<-as.character(complex[i])
  filename<-paste(folder,"phesant_",name,".txt",sep="")
  tmp<-read_tsv(file=filename) %>% select(-(isTraitOfInterest:Path),-(lower:upper))%>% filter(beta> -2)
  names(tmp)[4:5]<-paste0(name,c("_beta","_pvalue"))
  phesant<-full_join(phesant,tmp,by=c("varName","varType","n","resType","description"))
}
rm(tmp)
##split counts
phesant<-phesant %>% separate(n,c("control","case","total"))
## remove rows with high pvalue minimum
pvalue<-phesant %>% select(matches("pvalue")) %>% mutate(p_min=as.numeric(pmap(.,min,na.rm=TRUE))) %>% select(p_min)
phesant$p_min<-pvalue
phesant.trim<-phesant %>% filter(ALL_pvalue <1e-6,p_min>0) %>% mutate(across(control:total,~as.numeric(.)),case=ifelse(is.na(case),control,case))
#trim out some redundant codes
phesant.trim<-phesant.trim %>% filter(!str_detect(varName,"41204|41202|41200|41210|22190|22191"))
#select cancer and blood ICD codes
phesant.trim<-phesant.trim %>% filter(str_detect(varName,"41270#C|41270#D"))# | case>100000)
pvalue<-phesant.trim %>% filter(case>0)%>% select(description,matches("pvalue")) %>% mutate(across(ALL_pvalue:ComplexV_pvalue,~-log(.,10)))
rows.cor <- cor(t(pvalue %>% dplyr::select(-description,-ALL_pvalue)), use = "pairwise.complete.obs", method = "pearson")
rowlabels<-gsub("Diagnoses - ICD10: ","",pvalue$description)
colMain <- colorRampPalette(brewer.pal(8, "Reds"))(100)
pheatmap(( pvalue %>% dplyr::select(-description,-ALL_pvalue) %>% 
             mutate(across(DLOOP_pvalue:ComplexV_pvalue,~ifelse(.<2,0,.))) ),
         color=colMain,scale = "none",clustering_distance_row = as.dist(1 - rows.cor),cluster_cols = 0,
         labels_row =rowlabels,labels_col = gsub("Complex","Complex ",complex[-1]), display_numbers=TRUE,
         file = c("C:/Users/darking/Desktop/mt_phesant_heatmap.jpg" ),width = 9,height=4,angle_col = 45)