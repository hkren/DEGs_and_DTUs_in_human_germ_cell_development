###########
#libraries#
###########
install.packages("data.table")
if(!requireNamespace("BiocManager", quietly = T)){
  install.packages("BiocManager")
}
BiocManager::install(c("BiocParallel", "GenomicRanges", "Gviz", "rtracklayer", "stageR", "tximport"))

if(!requireNamespace("remotes", quietly = T)){
  install.packages("remotes")
}
remotes::install_github("TobiTekath/DTUrtle")

library(DTUrtle)

biocpar <- BiocParallel::MulticoreParam(32)

#############
#Preparation#
#############

#gtf Annotation to map transcripts to genes
tx2gene <- import_gtf(gtf_file = "gencode.v36.annotation.gtf")

#ensure that a one to one mapping between names and IDs exists
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)

tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

#move transcript and gene identifier columns to front
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

#import transcript-level quantification data, for example from Salmon
sampleMetadata <- fread("metadata.txt")

##############
#RNA samples #
##############
files <- Sys.glob(paste0("/salmon/366_", sampleMetadata$names, "_R1/quant.sf"))
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
cts <- import_counts(files = files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])

rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

#define grouping of samples
groups<-c("sco","sco","sco",
          "spg","spg","spg",
          "spz","spz","spz",
          "spd","spd","spd",
          "normal","normal","normal",
          "spg")
 
pd <- data.frame("id"=colnames(cts), "group"=groups, 
                 stringsAsFactors = F)

##############
#DTU Analysis#
##############

#dtu analysis for each pair of subsequent states
for(i in list(c("spz","spd"),c("spd","normal"),c("sco","normal"))){

  #use DRIMSeq for fitting a Dirichlet-multinomial model
  dturtle <- run_drimseq(counts = cts[!is.na(rownames(cts)),which(pd$group%in%i)], tx2gene = tx2gene, pd=pd, id_col = "id",
                         cond_col = "group", cond_levels = i, filtering_strategy = "bulk", 
                         BPPARAM = biocpar)
  
  save.image(file=paste("/workspace_dturtle", paste(i,collapse = " "),Sys.time(),".RData"))
  
  #run posthoc filtering and two-staged statistical correction with stageR
  dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)
  
  #highly flexible function to create a results data frame
  dturtle <- create_dtu_table(dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
                              add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))
  
  #change to results folder
  setwd("dturtle/")    
  
  #create plots, save them to disk and link them in the `dtu_table`.
  dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                     savepath = paste("images/",paste(i,collapse = " "),sep = ""), 
                                     add_to_table = "barplot",
                                     BPPARAM = biocpar)
  
  dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                   gtf = "gencode.v36.annotation.gtf", 
                                   genome = 'hg38', 
                                   one_to_one = T,
                                   savepath = paste("images/",paste(i,collapse = " "),sep = ""), 
                                   add_to_table = "transcript_view",
                                   BPPARAM = biocpar)
  
  #create interactive HTML-table from results data frame
  column_formatter_list <- list(
    "gene_qval" = table_pval_tile("white", "orange", digits = 3),
    "min_tx_qval" = table_pval_tile("white", "orange", digits = 3),
    "n_tx" = formattable::color_tile('white', "lightblue"),
    "n_sig_tx" = formattable::color_tile('white', "lightblue"),
    "max(Condition1-Condition2)" = table_percentage_bar('lightgreen', "#FF9999", digits=2))
  
  plot_dtu_table(dturtle = dturtle, savepath = paste("dturtle_",paste(i,collapse = " "),Sys.time(),".html"), 
                 column_formatters = column_formatter_list)
  
  #export proportiontable for DTUs
  proportionsTable<-get_proportion_matrix(dturtle$drim@counts[[dturtle$sig_gene[1]]])
  for (j in dturtle$sig_gene[2:length(dturtle$sig_gene)]){
    proportionsTable<-rbind(proportionsTable,get_proportion_matrix(dturtle$drim@counts[[j]]))
  }
  
  colnames(proportionsTable)<-paste(dturtle$meta_table_sample$condition,colnames(proportionsTable))
  
  proportionsTable<-data.frame("feature_id"=row.names(proportionsTable),proportionsTable)
  proportionsTable<-merge(y=proportionsTable, x=dturtle$meta_table_tx[,c("tx","transcript_id","transcript_type","gene","gene_type")], by.y="feature_id", by.x="tx",all.y=TRUE)
  
  tx_pvalue<-data.table(dturtle$FDR_table)
  names(tx_pvalue)[3:4]<-c("gene_FDR_value","transcript_FDR_value")
  
  proportionsTable<-merge(x=proportionsTable,y=tx_pvalue,by.x="tx",by.y="txID")
  proportionsTable$sig_transcript_dturtle<-proportionsTable$tx%in%dturtle$sig_tx
  
  fwrite(proportionsTable,paste("dturtle proportion table",
                                paste(i,collapse = " "),Sys.time(),".csv"))

}

#######################################
#proportion table: add effect strength#
#######################################
sco_spg_prop_table<-fread("dturtle proportion table sco spg.csv")
spg_spz_prop_table<-fread("dturtle proportion table spg spz.csv")
spz_spd_prop_table<-fread("dturtle proportion table spz spd.csv")
spd_normal_prop_table<-fread("dturtle proportion table spd normal.csv")

sco_spg_prop_table$'effect (mean(sco)-mean(spg))'<-rowMeans(sco_spg_prop_table[,list(sco.366_984_R1,sco.366_3425_R1,sco.366_653_R1)])-
  rowMeans(sco_spg_prop_table[,list(spg.366_3189_R1,spg.366_3309_R1,spg.366_1032_R1,spg.366_3445_R1)])
  
sco_spg_prop_table_filter_genelevel<-sco_spg_prop_table_filter[,list('max_positive_effect'=max(`effect (mean(sco)-mean(spg))`),
                                     'min_nagative_effect'=min(`effect (mean(sco)-mean(spg))`)),by=geneID]

sco_spg_prop_table_filter_genelevel[max_positive_effect<0,max_positive_effect:=NA]
sco_spg_prop_table_filter_genelevel[min_nagative_effect>0,min_nagative_effect:=NA]

sco_spg_prop_table_filter<-merge(sco_spg_prop_table_filter,sco_spg_prop_table_filter_genelevel,by="geneID")

fwrite(sco_spg_prop_table_filter,"dturtle proportion table filter sco spg.csv")

#
spg_spz_prop_table$'effect (mean(spg)-mean(spz))'<-rowMeans(spg_spz_prop_table[,list(spg.366_3189_R1,spg.366_3309_R1,spg.366_1032_R1,spg.366_3445_R1)])-
  rowMeans(spg_spz_prop_table[,list(spz.366_3103_R1,spz.366_529_R1,spz.366_658_R1)])

spg_spz_prop_table_filter_genelevel<-spg_spz_prop_table_filter[,list('max_positive_effect'=max(`effect (mean(spg)-mean(spz))`),
                                                                     'min_nagative_effect'=min(`effect (mean(spg)-mean(spz))`)),by=geneID]

spg_spz_prop_table_filter_genelevel[max_positive_effect<0,max_positive_effect:=NA]
spg_spz_prop_table_filter_genelevel[min_nagative_effect>0,min_nagative_effect:=NA]

spg_spz_prop_table_filter<-merge(spg_spz_prop_table_filter,spg_spz_prop_table_filter_genelevel,by="geneID")

fwrite(spg_spz_prop_table_filter,"dturtle proportion table filter spg spz.csv")

#
spz_spd_prop_table$'effect (mean(spz)-mean(spd))'<-rowMeans(spz_spd_prop_table[,list(spz.366_3103_R1,spz.366_529_R1,spz.366_658_R1)])-
  rowMeans(spz_spd_prop_table[,list(spd.366_932_R1,spd.366_3012_R1,spd.366_2038_R1)])

spz_spd_prop_table_filter_genelevel<-spz_spd_prop_table_filter[,list('max_positive_effect'=max(`effect (mean(spz)-mean(spd))`),
                                                                     'min_nagative_effect'=min(`effect (mean(spz)-mean(spd))`)),by=geneID]

spz_spd_prop_table_filter_genelevel[max_positive_effect<0,max_positive_effect:=NA]
spz_spd_prop_table_filter_genelevel[min_nagative_effect>0,min_nagative_effect:=NA]

spz_spd_prop_table_filter<-merge(spz_spd_prop_table_filter,spz_spd_prop_table_filter_genelevel,by="geneID")

fwrite(spz_spd_prop_table_filter,"dturtle proportion table filter spz spd.csv")

#
spd_normal_prop_table$'effect (mean(spd)-mean(normal))'<-rowMeans(spd_normal_prop_table[,list(spd.366_932_R1,spd.366_3012_R1,spd.366_2038_R1)])-
  rowMeans(spd_normal_prop_table[,list(normal.366_3178_R1,normal.366_927_R1,normal.366_1084_R1)])

spd_normal_prop_table_filter_genelevel<-spd_normal_prop_table_filter[,list('max_positive_effect'=max(`effect (mean(spd)-mean(normal))`),
                                                                     'min_nagative_effect'=min(`effect (mean(spd)-mean(normal))`)),by=geneID]

spd_normal_prop_table_filter_genelevel[max_positive_effect<0,max_positive_effect:=NA]
spd_normal_prop_table_filter_genelevel[min_nagative_effect>0,min_nagative_effect:=NA]

spd_normal_prop_table_filter<-merge(spd_normal_prop_table_filter,spd_normal_prop_table_filter_genelevel,by="geneID")

fwrite(spd_normal_prop_table_filter,"dturtle proportion table filter spd normal.csv")

#significant genes in all comparisons
sig_genes_in_all_comparisons<-merge(data.table("geneID"=sco_spg_prop_table_filter_genelevel$geneID,"sco_spg"=TRUE),
           data.table("geneID"=spg_spz_prop_table_filter_genelevel$geneID,"spg_spz"=TRUE),by="geneID",all=TRUE)

sig_genes_in_all_comparisons<-merge(sig_genes_in_all_comparisons,
                                    data.table("geneID"=spz_spd_prop_table_filter_genelevel$geneID,"spz_spd"=TRUE),by="geneID",all=TRUE)

sig_genes_in_all_comparisons<-merge(sig_genes_in_all_comparisons,
                                    data.table("geneID"=spd_normal_prop_table_filter_genelevel$geneID,"spd_normal"=TRUE),by="geneID",all=TRUE)

fwrite(sig_genes_in_all_comparisons,"signifikant genes in all comparisons 2021-04-23 .csv")

#################
#sankey diagramm#
#################
sco_spg_prop_table<-fread("dturtle proportion table sco spg.csv")
spg_spz_prop_table<-fread("dturtle proportion table spg spz.csv")
spz_spd_prop_table<-fread("dturtle proportion table spz spd.csv")
spd_normal_prop_table<-fread("dturtle proportion table spd normal.csv")

sco_spg_genes<-unique(sco_spg_prop_table$gene)
spg_spz_genes<-unique(spg_spz_prop_table$gene)
spz_spd_genes<-unique(spz_spd_prop_table$gene)
spd_normal_genes<-unique(spd_normal_prop_table$gene)

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
  value=c(sum(sco_spg_genes%in%spg_spz_genes),
          sum(sco_spg_genes%in%spg_spz_genes&
                sco_spg_genes%in%spz_spd_genes),
          sum(sco_spg_genes%in%spg_spz_genes&
                sco_spg_genes%in%spz_spd_genes&
                sco_spg_genes%in%spd_normal_genes),
          sum(sco_spg_genes%in%spg_spz_genes&
                !sco_spg_genes%in%spz_spd_genes&
                sco_spg_genes%in%spd_normal_genes),
          
          sum(!(sco_spg_genes%in%spg_spz_genes)&
                sco_spg_genes%in%spz_spd_genes),
          sum(!(sco_spg_genes%in%spg_spz_genes)&
                sco_spg_genes%in%spz_spd_genes&
                sco_spg_genes%in%spd_normal_genes),
          sum(!(sco_spg_genes%in%spg_spz_genes)&
                !(sco_spg_genes%in%spz_spd_genes)&
                sco_spg_genes%in%spd_normal_genes),
          
          sum(spg_spz_genes%in%spz_spd_genes&
                !spg_spz_genes%in%sco_spg_genes),
          sum(spg_spz_genes%in%spz_spd_genes&
                spg_spz_genes%in%spd_normal_genes&
                !spg_spz_genes%in%sco_spg_genes),
          
          sum(spg_spz_genes%in%spd_normal_genes&
                !(spg_spz_genes%in%spz_spd_genes)&
                !(spg_spz_genes%in%sco_spg_genes)),
          
          sum(spz_spd_genes%in%spd_normal_genes&
                !spz_spd_genes%in%sco_spg_genes&
                !spz_spd_genes%in%spg_spz_genes)
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

fwrite(links,"DTU_sankey.csv")

