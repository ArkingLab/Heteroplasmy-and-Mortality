#### Circos Plots ####
## This script to make circos plot for UKB heteroplasmy paper
## all files used to be saved in same folder on jhpce


#### Import files ####
library(data.table)
library(dplyr)
library(circlize)
library(R.utils)
library(ggplot2)
"%ni%" <- Negate("%in%")
library(bedr)

### circos 4 is the plot used for paper



#### circos 4 ####
## define everything in this chunk so it's neat


### make genes file ###
## downloaded from MitoHPC RefSeq files
dloop <- as.data.frame(fread('~/Desktop/UKB/circos_plot/DLOOP.bed'))
trna <- as.data.frame(fread('~/Desktop/UKB/circos_plot/TRN.bed'))
rrna <- as.data.frame(fread('~/Desktop/UKB/circos_plot/RNR.bed'))
cds <- as.data.frame(fread('~/Desktop/UKB/circos_plot/CDS.bed'))

dloop_mod <- data.frame(V1 = "chrM", V2 = c(1, 16023), V3 = c(576,16569), V4 = "Dloop")

genes <- rbind(dloop_mod, trna[,c("V1","V2","V3","V4")], rrna[,c("V1","V2","V3","V4")], cds)

genes <- genes[order(genes$V2),]

genes_fin <- subset(genes, V1 == "chrM")
genes_fin$color <- ifelse(genes_fin$V4 == "Dloop", "lightcoral", "lightgoldenrod")
genes_fin$color <- ifelse(genes_fin$V4 == "RNR1" | genes_fin$V4 == "RNR2", "thistle", genes_fin$color)
genes_fin$color <- ifelse(genes_fin$V4 == "ATP6" | genes_fin$V4 == "ATP8", "darkolivegreen3", genes_fin$color)
genes_fin$color <- ifelse(genes_fin$V4 == "CYTB", "darkseagreen2", genes_fin$color)
genes_fin$color <- ifelse(genes_fin$V4 == "COX1" | genes_fin$V4 == "COX2" | genes_fin$V4 == "COX3", "lightcyan2", genes_fin$color)
genes_fin$color <- ifelse(genes_fin$V4 == "ND1" | genes_fin$V4 == "ND2" | genes_fin$V4 == "ND3" | genes_fin$V4 == "ND4" | genes_fin$V4 == "ND4L" | genes_fin$V4 == "ND5" | genes_fin$V4 == "ND6", "lightblue", genes_fin$color)

# write.table(genes_fin, file="~/Desktop/UKB/circos_plot/genes.txt", row.names = F, quote = F, sep = "\t")


## redefined tracks below so everything in one place
df = data.frame(chr = "chrM", pos = seq(1,16569))
n = seq(1000, 16569, by = 1000)
labels = genes_fin$V4

#test <- as.data.frame(fread("~/Desktop/UKB/dna_nexus_results/per_pos_count_with_all_filter_0621.txt"))
test <- as.data.frame(fread("~/Desktop/UKB/per_pos_ref_alt_count_with_all_filter_0805.txt"))
test$Median_AF_Het[is.na(test$Median_AF_Het)] <- 0

## define syn, nonsyn, and stop in 1 col
test$coding <- ifelse(test$COMPLEX == "I" | test$COMPLEX == "III" | test$COMPLEX == "IV" | test$COMPLEX == "V", "coding", "noncoding")
test$col1 <- ifelse(test$coding == "coding" & test$mutation_nonsynonymous == "NONSYN", "NONSYN", NA)
test$col1 <- ifelse(test$coding == "coding" & is.na(test$mutation_nonsynonymous), "SYN", test$col1)
test$col1 <- ifelse(test$mutation_stop == "STOP" & !is.na(test$mutation_stop), "STOP", test$col1)

## define individual groups first
track3.1 <- subset(test, col1 == "SYN")
track3.1$chr <- "chrM"
track3.1$color <- ifelse(track3.1$COMPLEX == "I", "lightblue", "darkolivegreen3") ## complex I is lightblue, complex V is light green
track3.1$color <- ifelse(track3.1$COMPLEX == "III", "darkseagreen2", track3.1$color) ## cytb 
track3.1$color <- ifelse(track3.1$COMPLEX == "IV", "lightcyan2", track3.1$color)
track3.1$log_count_het <- log(track3.1$count_het + 1) ## because log(1) = 0

track3.2 <- subset(test, col1 == "NONSYN")
track3.2$chr <- "chrM"
track3.2$color <- ifelse(track3.2$COMPLEX == "I", "lightblue", "darkolivegreen3") ## complex I is lightblue, complex V is light green
track3.2$color <- ifelse(track3.2$COMPLEX == "III", "darkseagreen2", track3.2$color) ## cytb 
track3.2$color <- ifelse(track3.2$COMPLEX == "IV", "lightcyan2", track3.2$color)
track3.2$log_count_het <- log(track3.2$count_het + 1)

track3.3 <- subset(test, col1 == "STOP")
track3.3$chr <- "chrM"
track3.3$log_count_het <- log(track3.3$count_het + 1)
track3.3$color <- "red"

track3.4 <- subset(test, COMPLEX == "TRNA" | COMPLEX == "RRNA" | COMPLEX == "DLOOP")
track3.4$chr <- "chrM"
track3.4$color <- ifelse(track3.4$COMPLEX == "TRNA", "lightgoldenrod", "thistle")
track3.4$color <- ifelse(track3.4$COMPLEX == "DLOOP", "lightcoral", track3.4$color)
track3.4$log_count_het <- log(track3.4$count_het + 1)

## make the combined tracks
track3.5 <- rbind(track3.1[,c("chr","POS","log_count_het","color")],track3.4[,c("chr","POS","log_count_het","color")])
track3.6 <- rbind(track3.2[,c("chr","POS","log_count_het","color")],track3.3[,c("chr","POS","log_count_het","color")])

summary(track3.5$log_count_het)
summary(track3.6$log_count_het)

## make a dummy point to scale the track
scale_row <- data.frame(chr = "chrM", POS = c(1,2), log_count_het = c(8,0), color = "white")

meep <- rbind(track3.5, scale_row)
track3.5 <- meep

meep <- rbind(track3.6, scale_row)
track3.6 <- meep

## identify mutational deserts, two ways 
## first with 5 zeros in a row
d1 = df[which(df$pos %ni% test$POS),] ## get list of zeros
nrow(d1)

d2 <- data.frame(seqToIntervals(d1$pos))
d2$diff <- d2$to - d2$from
table(d2$diff)
# ggplot(d2, aes(x=diff)) +
#   geom_density() +
#   geom_rug() +
#   theme_classic()

## make 3 categories, 0,1,2+ which corespond to 1,2,3+ mutational desert sites in a row
d2$col1 = paste(d2$from,"-",d2$to,sep="")
zero = subset(d2, diff == 0)
one = subset(d2, diff == 1)
more = subset(d2, diff > 1)

l1 <- lapply(strsplit(zero$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)
l2 <- lapply(strsplit(one$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)
l3 <- lapply(strsplit(more$col1, "-"), function(x) Reduce(`:`, as.numeric(x))) %>% unlist(recursive = F)

track3.7 <- d1
track3.7$chr <- "chrM"
track3.7$color <- ifelse(track3.7$pos %in% l2, "gray57", "gray87") # everything in l1 should be yellow
track3.7$color <- ifelse(track3.7$pos %in% l3, "gray28", track3.7$color)
track3.7$col1 <- ifelse(track3.7$color == "gray28", 1, 0.33)
track3.7$col1 <- ifelse(track3.7$color == "gray57", 0.66, track3.7$col1)

scale_row2 <- data.frame(chr = "chrM", pos = 16565, col1 = 0, color = "white")

meep <- rbind(track3.7, scale_row2)
track3.7 <- meep


###### start plot ######
circos.clear()
circos.par("xaxis.clock.wise" = FALSE, "start.degree" = 90, "track.height" = 0.12)
circos.initialize(factors=df$chr, x=df$pos)

circos.track(factors=track3.5$chr, x = track3.5$POS, y = track3.5$log_count_het,
             panel.fun = function(x, y) {
               circos.lines(x, y, type='h',col = track3.5$color)
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 1, major.at = n)
             })

circos.track(factors=track3.6$chr, x = track3.6$POS, y = track3.6$log_count_het,
             panel.fun = function(x, y) {
               circos.lines(x, y, type='h',col = track3.6$color)
             })

circos.track(factors=track3.7$chr, x = track3.7$pos, y = track3.7$col1,
             panel.fun = function(x, y) {
               circos.lines(x, y, type='h',col = track3.7$color)
             })

circos.genomicTrack(genes_fin, stack = TRUE, genes_fin$V4, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, col = genes_fin$color, border = "black", ...)
                      #circos.genomicRect(region, value, col = genes_fin$color, border = "black", ytop = i + 0.3, ybottom = i - 0.3, ...)
                      #circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing = "inside", niceFacing = TRUE)
                      #circos.text(x,y,labels,cex = 0.8, facing = "inside")
                    },track.height=0.05,bg.border=F)
circos.genomicLabels(genes_fin, labels.column = 4, side = "inside", cex=0.70, connection_height = mm_h(2))





