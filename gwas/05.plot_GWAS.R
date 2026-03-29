library(GenomicRanges)
library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)
library(tidyverse)
library(ggbeeswarm)

#PCA
pca <- read_table2("Sw_merged_b1_b2_filtered_chr1.eigenvec", col_names = FALSE)
eigenval <- scan("Sw_merged_b1_b2_filtered_chr1.eigenval")
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
pca <- pca[-1,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
labs <- read.table("Sw_all_filt.fam")
pca$labs<-(labs[,6])
pca<-data.frame(pca)
pca$labs<-factor(pca$labs)
levels(pca$labs)<-c("Lw","Sw")
pops.pc1.chr1<-ggplot(pca, aes(x=as.numeric(PC1), y=labs,fill=labs,colour=labs))+
  scale_fill_manual(values=c('orange','purple'))+
  scale_colour_manual(values=c('orange','purple'))+
  geom_quasirandom(shape=21,size=1.5,aes(fill=labs),colour='#333',alpha=0.5)+theme_classic()+
  stat_summary(colour='black',shape=21,size=0.5,aes(fill=labs))+
  theme(panel.grid=element_blank(),axis.title.y=element_blank(),
        legend.position='none')+xlab(paste0("PC1 (",round(pve[1,2],digits=2),"%)"))

#fisher test
pca$inv<-FALSE
pca[pca$PC1>0,]$inv<-TRUE
fisher.test(table(pca$inv,pca$labs))

#GWAS
#first remove variants with P >= 0.1
Sw<-read.table("Sw_merged_b1_b2_filtered.assoc.txt",h=T)
head(Sw[order(Sw$p_lrt),])

Sw$padj<-p.adjust(Sw$p_lrt)
Sw<-Sw[Sw$p_lrt<0.1,]
write.csv(file="Sw_merged_b1_b2_filteredP0.1.csv",Sw)

#read in annotation info
tocgff<-read.table("~/Documents/StA/Sw_WGS/TOC.asm.scaffold.gene.gff3") %>% 
  mutate(V1=as.integer(gsub("scaffold_","",V1)),
         gene=gsub("ID=.*;.*Name=","",V9)) %>%
  filter(V3=="gene")

#read functional annotations
gene_annos<-read.table("~/Documents/StA/Sw_WGS/TOC.asm.scaffold.gene.SWISSPROT.blastp.top_hit.txt")
colnames(gene_annos)<-c("gene","score","anno")
gene_annos$gene<-gsub("\\.t.*","",gene_annos$gene)#remove transcript identifier

#keep top annotation per gene
gene_annos<-gene_annos[order(gene_annos$score,decreasing=TRUE),] %>% filter(!duplicated(gene_annos))

#remove accession info and spp ID from annotations for readability, then merge genome and functional annotation data frames
gene_annos$anno<-gsub(".*\\|","",gene_annos$anno)
tocgff.anno<-merge(tocgff,gene_annos,by='gene')
tocgff.anno$anno1<-gsub("_.*","",tocgff.anno$anno)

#create granges object used to identify overlapping regions (i.e. mapping variants to nearby genes)
genes<-GRanges(seqnames=tocgff.anno$V1,#chr/scaffolds
               ranges=IRanges(start=tocgff.anno$V4, end=tocgff.anno$V5))#start and end of gene feature
mcols(genes)$gene <- tocgff.anno$gene#add gene ID (e.g., TOC.11T.00352)
mcols(genes)$anno <- tocgff.anno$anno#add annotation (eg., UBCD1)

#read variants with P < 0.1

Sw<-read.csv("Sw_merged_b1_b2_filteredP0.1.csv",h=T)
Sw$chr<-as.integer(gsub("scaffold_","",Sw$chr))
Sw$sig<-Sw$padj<0.05
Sw<-annotation_function(Sw)

Sw.plot<-ggplot(Sw,aes(x=ps,y=(-log10(p_lrt)),colour=sig))+
  geom_rect(data=data.frame(chr=1),ymin=0,ymax=150,xmin=96e+06,xmax=246e+06,inherit.aes=FALSE,colour='#fff',fill='#f0f5ed',linewidth=0.05)+
  theme_classic()+facet_grid(.~(chr),space='free',scales='free')+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),legend.position='none',
        panel.spacing.x=unit(0.05, "lines"),
        strip.background=element_rect(colour='white',fill='white'),axis.line=element_line(),
        axis.title.x=element_blank(),plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))+
  geom_point(size=0.5)+scale_y_continuous(limits=c(1,10),expand=c(0,0))+
  scale_colour_manual(values=c("#ccc","#db6767"))+
  ylab("-log10(P)")

#create annotated df for sig variants in Sw-assoc region
tocgff.anno.sw<-tocgff.anno[tocgff.anno$V1==1 & tocgff.anno$V4>132436658 & tocgff.anno$V4<194378945,]
tocgff.anno.sw<-tocgff.anno.sw[!duplicated(tocgff.anno.sw$anno1),]
Sw.sig<-Sw[Sw$padj<0.05,]
window=50000
Sw.sig$window.start<-Sw.sig$ps-window
Sw.sig$window.end<-Sw.sig$ps+window
swsig<-GRanges(seqnames=Sw.sig$chr,ranges=IRanges(start=Sw.sig$window.start,end=Sw.sig$window.end))
swannos<-GRanges(seqnames=tocgff.anno.sw$V1,ranges=IRanges(start=tocgff.anno.sw$V4,end=tocgff.anno.sw$V5))
mcols(swannos)$anno<-tocgff.anno.sw$anno1
variants_within<-findOverlaps(swsig, swannos)
sw.anno<-data.frame(cbind(Sw.sig[variants_within@from,],tocgff.anno.sw[variants_within@to,]))
sw.anno<-sw.anno[!duplicated(sw.anno$anno1),]

Sw.plot.zoomed<-ggplot(Sw[Sw$chr==1 & Sw$ps>96e+06 & Sw$ps<246e+06,],aes(x=ps/1e+06,y=(-log10(p_lrt)),colour=sig))+
  geom_rect(data=data.frame(chr=1),ymin=0,ymax=150,xmin=96,xmax=246,inherit.aes=FALSE,colour='#fff',fill='#f0f5ed',linewidth=0.05)+
  theme_classic()+xlab("pos (Mb)")+
  theme(panel.grid=element_blank(),
        legend.position='none',
        panel.spacing.x=unit(0.05, "lines"),
        strip.background=element_rect(colour='white',fill='white'),axis.line=element_line())+
  geom_segment(data=test,aes(x=V4/1e+06,xend=V4/1e+06,y=1,yend=11),linetype='solid',colour="#db6767",alpha=1,linewidth=0.25)+
  geom_text_repel(data=sw.anno,
                  aes(label=anno1.1,x=V4/1e+06,y=11),
                  max.overlaps = 10,colour="#db6767",size=2,force_pull = 0.5,fontface = 'italic',
                  min.segment.length = 0,segment.size = 0.25,nudge_y = 1,nudge_x=0,segment.alpha=1)+
  scale_colour_manual(values=c("#ccc","#db6767"))+
  geom_point(size=1)+scale_y_continuous(limits=c(1,13.5),expand=c(0,0))+
  ylab("-log10(P)")



#
#R version 4.4.2 (2024-10-31)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sonoma 14.3

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: America/New_York
#tzcode source: internal

#attached base packages:
#[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] ggridges_0.5.6       ggbeeswarm_0.7.2     lubridate_1.9.4      forcats_1.0.0        purrr_1.0.4          readr_2.1.5         
# [7] tibble_3.2.1         tidyverse_2.0.0      tidyr_1.3.1          patchwork_1.3.0      ggplot2_4.0.2        dplyr_1.1.4         
#[13] stringr_1.5.1        GenomicRanges_1.56.2 GenomeInfoDb_1.40.1  IRanges_2.38.1       S4Vectors_0.42.1     BiocGenerics_0.50.0 

#loaded via a namespace (and not attached):
# [1] gtable_0.3.6            jsonlite_1.8.9          compiler_4.4.2          tidyselect_1.2.1        scales_1.4.0           
# [6] R6_2.5.1                XVector_0.44.0          generics_0.1.3          GenomeInfoDbData_1.2.12 tzdb_0.4.0             
#[11] pillar_1.10.1           RColorBrewer_1.1-3      rlang_1.1.5             stringi_1.8.4           S7_0.2.1               
#[16] timechange_0.3.0        cli_3.6.3               withr_3.0.2             magrittr_2.0.3          zlibbioc_1.50.0        
#[21] grid_4.4.2              rstudioapi_0.17.1       hms_1.1.3               beeswarm_0.4.0          lifecycle_1.0.4        
#[26] vipor_0.4.7             vctrs_0.6.5             glue_1.8.0              farver_2.1.2            httr_1.4.7             
#[31] tools_4.4.2             pkgconfig_2.0.3         UCSC.utils_1.0.0       
