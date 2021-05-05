###########
#libraries#
###########
library(limma)
library(BiocManager)
library(DESeq2)
library(edgeR)
library(data.table)
library(ggplot2)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(tximport)
library(knitr)
library("BiocParallel")
library(scales)
library(dendextend)
library(networkD3)
library(dplyr)

#############
#Preparation#
#############
palette(hue_pal()(5))
register(MulticoreParam(32))

#read in meta data
#table columns: run,	barcode,	names,	type,	state,	RIN
sampleMetadata <- fread("metadata.txt")
sampleMetadata$state<-factor(sampleMetadata$state,level=c("SCO","SPG","SPZ","SPD","normal"))

#create table of genes and ensembleIDs
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)

##############
#RNA samples #
##############

#prepare metadata
deseqSampleTable <- sampleMetadata

deseqSampleTable$files <- paste0("salmon/366_", deseqSampleTable$names, "_R1/quant.sf")
rownames(deseqSampleTable) <- deseqSampleTable$names

#import salmon files
txi <- tximport(deseqSampleTable$files, type ="salmon", tx2gene = t2g, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = deseqSampleTable,
                                   design = ~ state)

#keep transcripts with at least 10 counts over all samples
keep <- rowSums(counts(ddsTxi)) > 10
ddsTxi <- ddsTxi[keep,]

#create DEseq object
dds <- DESeq(ddsTxi)
colnames(dds)<-sampleMetadata$names

#pca
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("state"))
png(paste("PCA",Sys.Date(),".png",sep=""),
    units="in", width=5, height=9, res=300)
plotPCA(vsd, intgroup=c("state"))
dev.off()

#differential expression analysis for each pair of subsequent states
for(i in list(c("SCO","SPG"),c("SPG","SPZ"),c("SPZ","SPD"),c("SPD","normal"),c("SCO","normal"),c("SPG","normal"),
              c("SPZ","normal"))){

  #compute DEGs
  results_filter <- results(dds, contrast=c("state",i[1],i[2]), alpha=0.05)

  #remove NAs from results
  results_filter <- results_filter[!is.na(results_filter$padj),]
  
  #remove results with p$adj >= 0.05 from results
  results_filter_padj <- results_filter[results_filter$padj<0.05,]

  write.csv(results_filter_padj, 
            file=paste("DEGs_",paste(i,collapse="_"),sys.time(),".csv"))
  
  #extract counts of DEG genes for down stream analysis
  counts_normalized<-counts(dds, normalized=T)
  counts_raw<-counts(dds, normalized=F)

  counts_normalized<-counts_normalized[row.names(counts_normalized)%in%row.names(results_filter_padj),
                                     which(sampleMetadata$state%in%i)]
  colnames(counts_normalized)<-paste("normalized_counts",colnames(counts_normalized))
  
  counts_raw<-counts_raw[row.names(counts_raw)%in%row.names(results_filter_padj),
                       which(sampleMetadata$state%in%i)]
  colnames(counts_raw)<-paste("raw_counts",colnames(counts_raw))

  counts_all<-data.table("ensembl_gene_id"=row.names(counts_raw),data.table(counts_raw),data.table(counts_normalized),
                         as.data.table(results_filter_padj))

  t2g_table<-as.data.table(t2g[,2:3])
  t2g_table<-unique(t2g_table)
  counts_all<-t2g_table[counts_all, on = c(ensembl_gene_id='ensembl_gene_id')]
  
  write.csv(counts_all, 
          file=paste("Counts_DEGs_",paste(i,collapse="_"),sys.time(),".csv"))
}

##################################
#Sanky Diagram                   #
##################################

#read in DEGS
data_sco_spg<-fread("DEG_ SCO_SPG _counts_degs 2021-01-15 08:37:36 .csv")
data_spg_spz<-fread("DEG_ SPG_SPZ _counts_degs 2021-01-15 08:37:38 .csv")
data_spz_spd<-fread("DEG_ SPZ_SPD _counts_degs 2021-01-15 08:37:40 .csv")
data_spd_normal<-fread("DEG_ SPD_normal _counts_degs 2021-01-15 08:37:41 .csv")

#filter for log2Fold Change > 1 or < -1
data_sco_spg<-data_sco_spg[log2FoldChange<(-1)|log2FoldChange>1]
data_spg_spz<-data_spg_spz[log2FoldChange<(-1)|log2FoldChange>1]
data_spz_spd<-data_spz_spd[log2FoldChange<(-1)|log2FoldChange>1]
data_spd_normal<-data_spd_normal[log2FoldChange<(-1)|log2FoldChange>1]

links <- data.frame(
  source=c("sco vs spg","spg vs spz","spz vs spd","spg vs spz",
           "sco vs spg","spz vs spd",
           "sco vs spg",
           "spg vs spz","spz vs spd",
           "spg vs spz",
           "spz vs spd"),
  target=c("spg vs spz","spz vs spd","spd vs normal","spd vs normal",
           "spz vs spd","spd vs normal",
           "spd vs normal",
           "spz vs spd","spd vs normal",
           "spd vs normal",
           "spd vs normal"), 
  group=as.character(c(1,1,1,1,2,2,3,4,4,5,6)),
  value=c(sum(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name),
          sum(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name),
          sum(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name),
          sum(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                !data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name),
          
          sum(!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name),
          sum(!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name),
          sum(!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                !(data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name),
          
          sum(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name),
          sum(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name),
          
          sum(data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                !(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name)&
                !(data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name)),

          sum(data_spz_spd$external_gene_name%in%data_spd_normal$external_gene_name&
                !data_spz_spd$external_gene_name%in%data_sco_spg$external_gene_name&
                !data_spz_spd$external_gene_name%in%data_spg_spz$external_gene_name)
  ),
  sourc_up_regulation=c(data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&log2FoldChange>0,.N],
               data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&log2FoldChange>0,.N],
               data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange>0,.N],
               data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                !data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange>0,.N],
          
               data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&log2FoldChange>0,.N],
               data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange>0,.N],
               data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                !(data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name)&
                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange>0,.N],
          
          data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name&log2FoldChange>0,.N],
          data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name&log2FoldChange>0,.N],
          
          data_spg_spz[data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                !(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name)&
                !(data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name)&log2FoldChange>0,.N],
          
          data_spz_spd[data_spz_spd$external_gene_name%in%data_spd_normal$external_gene_name&
                !data_spz_spd$external_gene_name%in%data_sco_spg$external_gene_name&
                !data_spz_spd$external_gene_name%in%data_spg_spz$external_gene_name&log2FoldChange>0,.N]
  ),
  source_down_regulation=c(data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&log2FoldChange<0,.N],
                  data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                 data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&log2FoldChange<0,.N],
                  data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                 data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                 data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange<0,.N],
                  data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                 !data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                 data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange<0,.N],
                  
                  data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                 data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&log2FoldChange<0,.N],
                  data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                 data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                 data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange<0,.N],
                  data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                 !(data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name)&
                                 data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name&log2FoldChange<0,.N],
                  
                  data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                 !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name&log2FoldChange<0,.N],
                  data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                 data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                 !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name&log2FoldChange<0,.N],
                  
                  data_spg_spz[data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                 !(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name)&
                                 !(data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name)&log2FoldChange<0,.N],
                  
                  data_spz_spd[data_spz_spd$external_gene_name%in%data_spd_normal$external_gene_name&
                                 !data_spz_spd$external_gene_name%in%data_sco_spg$external_gene_name&
                                 !data_spz_spd$external_gene_name%in%data_spg_spz$external_gene_name&log2FoldChange<0,.N]
  ),
  target_up_regulation=c(data_spg_spz[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spz_spd[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                          data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                          data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                          data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                          !data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                          data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           
                           data_spz_spd[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                          data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                          data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                          data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                          !(data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name)&
                                          data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           
                           data_spz_spd[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                          !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           data_spd_normal[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                          data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                          !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name,external_gene_name]&log2FoldChange>0,.N],
                           
                           data_spd_normal[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                          !(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name)&
                                          !(data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name),external_gene_name]&log2FoldChange>0,.N],
                           
                           data_spd_normal[external_gene_name%in%data_spz_spd[data_spz_spd$external_gene_name%in%data_spd_normal$external_gene_name&
                                          !data_spz_spd$external_gene_name%in%data_sco_spg$external_gene_name&
                                          !data_spz_spd$external_gene_name%in%data_spg_spz$external_gene_name,external_gene_name]&log2FoldChange>0,.N]
  ),
  target_down_regulation=c(data_spg_spz[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spz_spd[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                                                             data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                                                                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                                                                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name&
                                                                                !data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                                                                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           
                           data_spz_spd[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                                                             data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                                                                data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name&
                                                                                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spd_normal[external_gene_name%in%data_sco_spg[!(data_sco_spg$external_gene_name%in%data_spg_spz$external_gene_name)&
                                                                                !(data_sco_spg$external_gene_name%in%data_spz_spd$external_gene_name)&
                                                                                data_sco_spg$external_gene_name%in%data_spd_normal$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           
                           data_spz_spd[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                                                             !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           data_spd_normal[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name&
                                                                                data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                                                                !data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name,external_gene_name]&log2FoldChange<0,.N],
                           
                           data_spd_normal[external_gene_name%in%data_spg_spz[data_spg_spz$external_gene_name%in%data_spd_normal$external_gene_name&
                                                                                !(data_spg_spz$external_gene_name%in%data_spz_spd$external_gene_name)&
                                                                                !(data_spg_spz$external_gene_name%in%data_sco_spg$external_gene_name),external_gene_name]&log2FoldChange<0,.N],
                           
                           data_spd_normal[external_gene_name%in%data_spz_spd[data_spz_spd$external_gene_name%in%data_spd_normal$external_gene_name&
                                                                                !data_spz_spd$external_gene_name%in%data_sco_spg$external_gene_name&
                                                                                !data_spz_spd$external_gene_name%in%data_spg_spz$external_gene_name,external_gene_name]&log2FoldChange<0,.N]
  )
)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, LinkGroup = "group")
p

png(paste("sankey_diagram_all",Sys.Date(),".png",sep=""),
    units="in", width=5, height=9, res=300)
p
dev.off()

fwrite(links,"DEGs_sanky.csv")

#common DEGs in all comparisons
sig_genes_in_all_comparisons<-merge(data.table("externalGeneName"=data_sco_spg$external_gene_name,"sco_spg"=TRUE),
                                    data.table("externalGeneName"=data_spg_spz$external_gene_name,"spg_spz"=TRUE),by="externalGeneName",all=TRUE)

sig_genes_in_all_comparisons<-merge(sig_genes_in_all_comparisons,
                                    data.table("externalGeneName"=data_spz_spd$external_gene_name,"spg_spz"=TRUE),by="externalGeneName",all=TRUE)

sig_genes_in_all_comparisons<-merge(sig_genes_in_all_comparisons,
                                    data.table("externalGeneName"=data_spd_normal$external_gene_name,"spg_spz"=TRUE),by="externalGeneName",all=TRUE)

fwrite(sig_genes_in_all_comparisons,"sankey signifikant genes in all comparisons.csv")
