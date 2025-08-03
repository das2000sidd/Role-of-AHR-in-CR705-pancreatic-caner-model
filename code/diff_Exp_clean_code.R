setwd("~/Desktop/PhD_Project_related/CR705_tumor_RNA_seq/Counts file")


library("DESeq2");
library("edgeR");
library("dplyr");
library("ggplot2");
library("pheatmap");
library("RColorBrewer");
library("Glimma");
library("org.Mm.eg.db");

## Reading in gene expression counts table for WT
wt_samples <- setNames(lapply(
	c(paste0("G3M",c(1,3,4,6,7,9))), 
	function(x) read.table(paste0(x, ".htseq.counts.txt"), 
	header = FALSE, 
	row.names = 1)), 
	c(paste0("G3M",
	c(1,3,4,6,7,9))
	)
);

combined_wt <- do.call(cbind, wt_samples);

## Set colnames of counts to WT1 to WT5
colnames(combined_wt) <- paste0("WT",1:6);

## Reading in gene expression counts table for mutant
mutant_samples <- setNames(lapply(
	c(paste0("G4M",c(3,4,5,7,9))), 
	function(x) read.table(paste0(x, ".htseq.counts.txt"), 
	header = FALSE, 
	row.names = 1)), 
	c(paste0("G4M",c(3,4,6,7,9))
	)
);


combined_ahrko <- do.call(cbind, 
	mutant_samples
	);

## Set colnames of counts to AHRKO1 to AHRKO5
colnames(combined_ahrko) <- paste0("AHRKO",1:5);


## Combine WT and AHRKO count in single table
all_counts_combined <- cbind(combined_wt,
	combined_ahrko[rownames(combined_wt),]
    );


## Define the grouping variable for sample in order of counts
condition_df <- as.data.frame(c(rep("WT",6),
rep("AHRKO",5)
));

## Set grouping variable as a factor
colnames(condition_df) <- "condition";
condition_df$condition <- as.factor(condition_df$condition);



## Define DESeq2 object for differential expression
dds <- DESeqDataSetFromMatrix(countData = all_counts_combined,
	colData = condition_df,
	design = ~ condition
	);



## Estimate sample size factor to account for difference in library size and composition
dds <- estimateSizeFactors(dds);

## Generate the variance stabilisation transformed object for plotting PCA
vsd <- vst(dds, 
	blind = FALSE
	);

head(assay(vsd), 3);

###PCA plot
pca.plot <- plotPCA(vsd, 
	intgroup = c("condition")) + ggtitle("PCA by genotype");


## Save PCA plot to file
pdf(file="CR705_MDS_plot.pdf",
	width = 10,
	height = 10);

pca.plot;

dev.off();



## Fit the linear model using DESeq command upon which contrasts will be run eventually
dds <- DESeq(dds);


## Defining contrast and pulling out the differential expresison table suing results
## AHRKO is comparison group and WT is reference group
res_AHRKO_vs_WT <- results(dds, 
	contrast = c("condition","AHRKO","WT"),
	pAdjustMethod = "BH",
	cooksCutoff = FALSE,
	independentFiltering = FALSE
	);


res_AHRKO_vs_WT_df <- as.data.frame(res_AHRKO_vs_WT);

res_AHRKO_vs_WT_df$Ensembl <- rownames(res_AHRKO_vs_WT_df);

res_AHRKO_vs_WT_df <- res_AHRKO_vs_WT_df[! is.na(res_AHRKO_vs_WT_df$padj),];

## Annoating ensemble IDs using org.Mm.eg.db package to gene gene symbol, name via entrez
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
                                      column="GENENAME"
                                      );


res_AHRKO_vs_WT_df$Entrez <- as.character(res_AHRKO_vs_WT_df$Entrez);

res_AHRKO_vs_WT_df$Symbol <- as.character(res_AHRKO_vs_WT_df$Symbol);

res_AHRKO_vs_WT_df$Genename <- as.character(res_AHRKO_vs_WT_df$Genename);


## Pulling out significantly up and down genes
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

## Write differential expression table to file for plotting
write.csv(res_AHRKO_vs_WT_df,
          file="AHRKO_vs_WT_differential_Expression_Table.csv",
          col.names = T,
          row.names = F,
          quote = F
          );

## The normalised counts which will be used for plotting MA plot and Z score ehatmap
normalised_counts=counts(dds, 
                         normalized=T 
                         );

write.csv(normalised_counts,
          file="AHRKO_vs_WT_Normalised_Expression_Table.csv",
          col.names = T,
          row.names = T,
          quote = F
          )


## save library info
sink("Differential_expression_sesion_info.txt");
sessionInfo();
sink();



