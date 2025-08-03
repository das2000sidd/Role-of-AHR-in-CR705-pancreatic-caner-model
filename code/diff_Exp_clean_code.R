setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")


library("DESeq2");
library("edgeR");
library("dplyr");
library("ggplot2");
library("pheatmap");
library("RColorBrewer");
library("Glimma");
library("org.Mm.eg.db");


wt_samples <- setNames(lapply(
                    			c(paste0("G3M",c(1,3,4,6,7,9))), 
                    			function(x) read.table(paste0(x, ".htseq.counts.txt"), 
                    								   header = FALSE, 
                    								   row.names = 1
                    								   )
                    				), 
                    			c(paste0("G3M",
                    					c(1,3,4,6,7,9)
                    					)
                    			 )
                    	);

combined_wt <- do.call(cbind, wt_samples);

colnames(combined_wt) <- paste0("WT",1:6);

mutant_samples <- setNames(lapply(
                       			   c(paste0("G4M",c(3,4,5,7,9))), 
                       			   function(x) read.table(paste0(x, ".htseq.counts.txt"), 
                       									 header = FALSE, 
                       									 row.names = 1
                       									 )
                       						), 
                      				c(paste0("G4M",
                      						c(3,4,6,7,9)
                      						)
                      			      )
                      		);


combined_ahrko <- do.call(cbind, 
                          mutant_samples
                          );

colnames(combined_ahrko) <- paste0("AHRKO",1:5);



all_counts_combined <- cbind(combined_wt,
                          	 combined_ahrko[rownames(combined_wt),]
                          	)



condition_df <- as.data.frame(c(rep("WT",6),
                             	rep("AHRKO",5)
                             	)
                           	);


colnames(condition_df) <- "condition";

condition_df$condition <- as.factor(condition_df$condition);



dds <- DESeqDataSetFromMatrix(countData = all_counts_combined,
                             colData = condition_df,
                             design = ~ condition
                             );




dds <- estimateSizeFactors(dds);


vsd <- vst(dds, 
           blind = FALSE
           );

head(assay(vsd), 3);

###PCA plot***

pca.plot <- plotPCA(vsd, 
        intgroup = c("condition")
        ) + 
        ggtitle("PCA by genotype");


## MDS plot using the VST data
pdf(file="CR705_MDS_plot.pdf",width = 10,height=10);
pca.plot;
dev.off();



## Linear mdoel fitting
dds <- DESeq(dds);


## Calling result to get differential expression table
res_AHRKO_vs_WT <- results(dds, 
                           contrast = c("condition",
                                        "AHRKO",
                                        "WT"
                                        ),
                           pAdjustMethod = "BH",
                           cooksCutoff = FALSE,
                           independentFiltering = FALSE
                           );


res_AHRKO_vs_WT_df <- as.data.frame(res_AHRKO_vs_WT);

res_AHRKO_vs_WT_df$Ensembl <- rownames(res_AHRKO_vs_WT_df);

res_AHRKO_vs_WT_df <- res_AHRKO_vs_WT_df[! is.na(res_AHRKO_vs_WT_df$padj),];


res_AHRKO_vs_WT_df$Entrez <- mapIds(org.Mm.eg.db, 
                                    res_AHRKO_vs_WT_df$Ensembl,
                                    keytype="ENSEMBL", 
                                    column="ENTREZID",
                                    multiVals = "first"
                                  );

res_AHRKO_vs_WT_df$Symbol <- mapIds(org.Mm.eg.db, 
                                    res_AHRKO_vs_WT_df$Entrez,
                                    keytype="ENTREZID", 
                                    column="SYMBOL",
                                    multiVals = "first"
                                    );


res_AHRKO_vs_WT_df$Genename <- mapIds(org.Mm.eg.db, 
                                      res_AHRKO_vs_WT_df$Entrez,
                                      keytype="ENTREZID", 
                                      column="GENENAME");


res_AHRKO_vs_WT_df$Entrez <- as.character(res_AHRKO_vs_WT_df$Entrez);

res_AHRKO_vs_WT_df$Symbol <- as.character(res_AHRKO_vs_WT_df$Symbol);

res_AHRKO_vs_WT_df$Genename <- as.character(res_AHRKO_vs_WT_df$Genename);



up_keeping_AHR_KO_1 <- subset(res_AHRKO_vs_WT_df,
                              res_AHRKO_vs_WT_df$log2FoldChange > 1 & 
                                res_AHRKO_vs_WT_df$padj < 0.01
                              );

dn_keeping_AHR_KO_1 <- subset(res_AHRKO_vs_WT_df,
                              res_AHRKO_vs_WT_df$log2FoldChange < -1 & 
                                res_AHRKO_vs_WT_df$padj < 0.01
                              );


print("Number of up genes");
nrow(up_keeping_AHR_KO_1);

print("Number of down genes");
nrow(dn_keeping_AHR_KO_1)

write.csv(res_AHRKO_vs_WT_df,
          file="AHRKO_vs_WT_differential_Expression_Table.csv",
          col.names = T,
          row.names = F,
          quote = F
          );


normalised_counts=counts(dds, 
                         normalized=T 
                         );

write.csv(normalised_counts,
          file="AHRKO_vs_WT_Normalised_Expression_Table.csv",
          col.names = T,
          row.names = T,
          quote = F
          )



sink("Differential_expression_sesion_info.txt");
sessionInfo();
sink();



