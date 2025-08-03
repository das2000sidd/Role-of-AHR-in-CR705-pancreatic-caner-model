setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")


library(dplyr);
library(clusterProfiler);
library(msigdbr);
library(org.Mm.eg.db);
library(magrittr);
library(org.Mm.eg.db);



`%ni%` <- Negate(`%in%`);

## Read in differential expression table
wt_vs_ahrko <- read.csv(file="WT_vs_AHRKO_All_genes_differential_Expression_Table.csv",
                        header = T,
                        stringsAsFactors = FALSE
                        );

wt_vs_ahrko$Entrez <- as.character(wt_vs_ahrko$Entrez);


## Get the gene names
wt_vs_ahrko$Genename <- mapIds(org.Mm.eg.db,
                               wt_vs_ahrko$Entrez,
                               keytype="ENTREZID", 
                               column="GENENAME"
                               );

wt_vs_ahrko$Genename <- as.character(wt_vs_ahrko$Genename);

## Remove riken, predicted and pseudogene before enrichment
ahrko_noGm_riken <- wt_vs_ahrko[- grep("RIKEN",wt_vs_ahrko$Genename),];

ahrko_noGm_riken <- ahrko_noGm_riken[- grep("predicted",ahrko_noGm_riken$Genename),];

ahrko_noGm_riken <- ahrko_noGm_riken[- grep("Riken",ahrko_noGm_riken$Genename),];

ahrko_noGm_riken <- ahrko_noGm_riken[- grep("pseudogene",ahrko_noGm_riken$Genename),];



## subset to significant genes only
wt_vs_ahrko_s=subset(ahrko_noGm_riken,
                     abs(ahrko_noGm_riken$log2FoldChange) > 1 & 
                     ahrko_noGm_riken$padj < 0.01
                     ); ## 1597



## pull out up and down genes separately
wt_vs_ahrko_up <- subset(wt_vs_ahrko_s,
                         wt_vs_ahrko_s$log2FoldChange > 0
                         );

wt_vs_ahrko_dn <- subset(wt_vs_ahrko_s,
                         wt_vs_ahrko_s$log2FoldChange < 0
                         );


write.csv(wt_vs_ahrko[,-c(10)],
          file="WT_vs_AHRKO_clean_list_genes.csv",
          col.names = T,
          row.names = F,
          quote = F
          );

write.csv(wt_vs_ahrko_up[,-c(10)],
          file="WT_vs_AHRKO_up_list_genes.csv",
          col.names = T,
          row.names = F,
          quote = F
          );

write.csv(wt_vs_ahrko_dn[,-c(10)],
          file="WT_vs_AHRKO_down_list_genes.csv",
          col.names = T,
          row.names = F,
          quote = F
          );



## Enrichment analysis
msigdbr_species()

mm_msigdb_df <- msigdbr(species = "Mus musculus");


#Filter the human data frame to the GO and Reactome pathway terms that are included in the
# curated gene sets
hs_GO_df <- mm_msigdb_df %>%
  dplyr::filter(
                gs_collection == "C5", 
                gs_subcollection %in% c("GO:BP","GO:CC","GO:MF") 
                );


hs_Reactome_df <- mm_msigdb_df %>%
  dplyr::filter(
                gs_collection == "C2", 
                gs_subcollection %in% c("CP:REACTOME") 
               );



## Declare background genes for enrichment
background_genes=c(wt_vs_ahrko$Ensembl);

background_genes=unique(background_genes);

## GO terms enrichment
GO_ora_results <- enricher(
                            gene = wt_vs_ahrko_s$Ensembl, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe = background_genes,
                            TERM2GENE = dplyr::select(
                                                      hs_GO_df,
                                                      gs_name,
                                                      ensembl_gene
                                                      )
                            );



View(GO_ora_results@result)


enrich_plot_go <- enrichplot::dotplot(GO_ora_results, 
                                   showCategory = 15,
                                   font.size = 6,
                                   title = "GO terms for genes up+down in AHRKO vs WT, CR705",
                                   orderBy = "p.adjust", 
                                   decreasing = FALSE
                                   );
                                   
enrich_plot_go;



tiff(file="AHRKO_vs_WT_CR705_all_sig_genes_GO_enrichment.tiff",
     res=300,
     height = 1500,
     width = 3000
     );

grid.draw(enrich_plot_go);

dev.off()


## Reactome terms enrichment
Reactome_ora_results <- enricher(
                                  gene = wt_vs_ahrko_s$Ensembl, 
                                  pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH",
                                  universe = background_genes,
                                  TERM2GENE = dplyr::select(
                                                            hs_Reactome_df,
                                                            gs_name,
                                                            ensembl_gene
                                                            )
                                  );



View(Reactome_ora_results@result)


enrich_plot_reac <- enrichplot::dotplot(Reactome_ora_results, 
                                   showCategory = 15,
                                   font.size = 6,
                                   title = "Reactome terms for genes up+down in AHRKO vs WT, CR705",
                                   orderBy = "p.adjust", 
                                   decreasing = FALSE
                                   );
                                   
enrich_plot_reac;



pdf(file = "AHRKO_vs_WT_CR705_significant_genes_Reactome_enrichment.pdf",
    width = 8,
    height = 8
    );

grid.draw(enrich_plot_reac);

dev.off();




sink("Enrichment_sesion_info.txt");
sessionInfo();
sink();






