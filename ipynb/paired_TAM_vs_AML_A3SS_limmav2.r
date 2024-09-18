# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .r
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.3
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# Paired TAM (transcient abnormal myleoproliferation) vs AML (acute myeloid leukemia) analysis of patients with the co-occuring condition of Down Syndrome using DESeq2 on IJC counts obtained from rMATS analysis.
#
# Using a matrix constructed from Kids First Workflow V4 done on single runs, a series of scripts were created and are stored in this repository.  For each of the splicing types, all the runs considered for analysis are pooled and normalized to have a non-redundant set of splicing events.  A matrix is then constructed for each of the samples to be analyzed.  
#
# Each splicing type has a bed file for visualizaiton in UCSC Genome browser of all the events, as well as created a matrix of the single runs normalized to the non-redundant union of files.  Both the source and the normalized bed file are available to ensure interprebility of results. 
#
# Using associative arrays in an awk script, it was a rapid way to transform the individual counts from each of the individual runs into a matrix that facilitated analysis.
#
# Using annotations obtained from the rMATS run that provided the coordinates of each of the splicing events as well as the gene that the junctions came from and the count of the reads that overlapped the junctions.   
#  
# Limma in this notebook is used to perform analysis of these junction counts provided by the rMATS routine.  Using these counts as junction expression.
#
# Between the splicing event differences and the expression differences, between paired samples, biological functional differences may be obtained.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")


# The viper tools are installed from the command line using `conda install` functions
#

library(mixtools)

library(bcellViper)

library(MASS)

library("viper")

BiocManager::install("dplyr")

library(Glimma)
library(dplyr)
library(edgeR)

setwd("/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/paired.TAM.AMLv2/A3SS_calculate")


getwd()


cts <- as.matrix(read.csv("A3SS.IJC.w.coordinates.matrix.csv",sep=",",row.names="ID"))

cts[1:3,11:dim(cts)[2]]

featureData <- data.frame(cts[,1:10])
featureData[1:3,]

featureData <- featureData[,c(1,2)]

head(featureData,2)

cts <- data.matrix(cts[,11:20])
mode(cts) <- "integer"
is.integer(cts)

dim(cts)
head(cts,2)

colnames(cts) <- c("PAUVKY.03A","PAUVKY.40A","PAWHSD.03A","PAWHSD.40A","PAWSNZ.03A","PAWSNZ.40A","PAUTLA.03A","PAUTLA.40A","PAVUDU.03A","PAVUDU.40A")

# The PAWHSD samples are not TAM and AML but in fact TAM and TAM - the resulting heatmap when included showed they clustered together.  We will eliminate them from subsequent analyses.

cts <- cts[,-c(3,4)]

head(cts,2)

# +
# Condition 1: Rows with count > 1000 in columns 1, 3, 5, 7
TAM_rows_condition <- rowSums(cts[, c(1, 3, 5, 7)] > 1000) > 0
TAM_matrix<- cts[TAM_rows_condition, c(1, 3, 5, 7)]

# Condition 2: Rows with count > 1000 in columns 2, 4, 6, 8
AML_rows_condition <- rowSums(cts[, c(2, 4, 6, 8)] > 1000) > 0
AML_matrix <- cts[AML_rows_condition, c(2, 4, 6, 8)]

# Combine the sub-matrices by keeping rows that satisfy either condition 1 or condition 2
final_matrix <- cts[AML_rows_condition | TAM_rows_condition, ]

# View the dimensions of the resulting matrices
dim(TAM_matrix)
dim(AML_matrix)
dim(final_matrix)

# -

head(TAM_matrix)
head(AML_matrix)
head(final_matrix)

# +
cts <- final_matrix

featureData <- featureData[rownames(cts),]
# -

dim(cts)
dim(featureData)

coldata <- read.csv("/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/paired.TAM.AMLv2/design_matrix.csv",row.names=1)

coldata


coldata <- coldata[,c("patient","condition")]
coldata$condition <- factor(coldata$condition)
coldata$patient <- factor(coldata$patient)

rownames(coldata)

rownames(coldata) <-sub("-",".",rownames(coldata))

colnames(cts)

all(rownames(coldata) %in% colnames(cts))

dim(cts)
head(cts,4)
mode(cts) <- "integer"
is.integer(cts)

# lets look at limma/voom
BiocManager::install("limma")

BiocManager::install("statmod")

library(limma)
library(edgeR)
library(statmod)

# making a counts matrix
dge <- DGEList(counts=cts)

colnames(dge)

head(dge,2)

# ## Explanation
# 1. **Library loading:** We load the limma and edgeR libraries required for the analysis.
# 2. **Group factor:** We create a factor group based on the condition column in coldata.
# 3. **Design matrix:** We create a design matrix using model.matrix where each column corresponds to a condition.
# 4. **Contrast matrix:** We define a contrast matrix using makeContrasts for the comparisons of interest.
# 5. **DGEList object:** We convert the count data to a DGEList object and normalize it.
# 6. **Voom transformation:** We apply the voom transformation to the normalized data.
# 7. **Linear model fitting:** We fit a linear model to the transformed data using lmFit.
# 8. **Contrast fitting:** We apply the contrasts to the fitted model using contrasts.fit.
# 9. **eBayes:** We compute the statistics using eBayes.
# 10. **TopTable:** We extract the top differentially expressed genes for each comparison using topTable.

# +
# Assuming cts is your count matrix and coldata is your sample information
# Make sure coldata$name matches the column names of cts
all(coldata$name %in% colnames(cts)) # This should return TRUE

# Create a factor for the conditions
group <- factor(coldata$condition)

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Print the design matrix to verify
print(design)

# Define the contrast matrix for the comparisons
# We want to compare TAM vs preAML, TAM vs AML, preAML vs AML
contrast_matrix <- makeContrasts(
  TAM_vs_preAML = TAM - preAML,
  TAM_vs_AML = TAM - AML,
  preAML_vs_AML = preAML - AML,
  levels = design
)

# Print the contrast matrix to verify
print(contrast_matrix)

# Convert the count data to a DGEList object
dge <- DGEList(counts = cts)

# Normalize the data using the TMM method
dge <- calcNormFactors(dge)

# Apply the voom transformation
v <- voom(dge, design)

# Fit the linear model
fit <- lmFit(v, design)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute the statistics
fit2 <- eBayes(fit2)

# -

design

# normalize and filter
keep          <-filterByExpr(dge, design)

is.logical(keep)
sum(keep==TRUE)

dge          <- dge         [keep,,keep.lib.size=FALSE]

# apply scale normalization
dge          <- calcNormFactors(dge)

# MDS Plot - can we separate the samples well?
logCPM <- cpm(dge, log=TRUE, prior.count=3)
plotMDS(logCPM,labels=coldata$condition,top=10, col=c(rep(c("red","black"),3)))

head(logCPM)

head(dge,2)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
de_results <- topTable(fit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(featureData[lookup,2])
head(featureData[lookup,2])

# +
# There are too many values - lets reduce the size a bit more
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 5
adjusted_pvalue_threshold <- 0.05

# Select genes that meet both fold change and adjusted p-value criteria
significant_genes <- de_results[
  abs(de_results$logFC) > fold_change_threshold &
  de_results$adj.P.Val < adjusted_pvalue_threshold,
]
dim(significant_genes)

# +
lookup <- rownames(significant_genes)
df <- as.data.frame(coldata[,c("condition","patient")])

significant_expression <- dge[lookup,]
dim(significant_expression)
length(lookup)

# +
library("pheatmap")

significant_out <- pheatmap(significant_expression, 
                            cluster_rows5=TRUE, 
                            show_rownames=FALSE,
                            cluster_cols=TRUE, 
                            annotation_col=df, 
                            scale="row",
                            clustering_method = "ward.D2",
                            clustering_distance_cols = "minkowski", 
                            clustering_distance_rows = "minkowski" )
# -

# weighting 
v <- voom(dge, plot=TRUE, normalize="quantile")

vfit <- lmFit(v, design)
vfit <- eBayes(vfit, trend=TRUE)
de_results <- topTable(vfit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(featureData[lookup,2])
head(featureData[lookup,2])

# +
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 9
adjusted_pvalue_threshold <- 0.05

# Select genes that meet both fold change and adjusted p-value criteria
significant_genes <- de_results[
  abs(de_results$logFC) > fold_change_threshold &
  de_results$adj.P.Val < adjusted_pvalue_threshold,
]
dim(significant_genes)
# -

lookup <- rownames(significant_genes)
significant_expression <- dge[lookup,]

significant_out <- pheatmap(significant_expression, 
                            cluster_rows5=TRUE, 
                            show_rownames=FALSE,
                            cluster_cols=TRUE, 
                            annotation_col=df, 
                            scale="row",
                            clustering_method = "ward.D2",
                            clustering_distance_cols = "minkowski", 
                            clustering_distance_rows = "minkowski" )

featureData[head(rownames(significant_expression),5),2]

top_gene_list <- as.matrix(featureData[rownames(significant_expression),2])
length(top_gene_list)

top_significant_genes <- dge[rownames(significant_genes),]

# +
start=1
stop=length(top_gene_list)
date="2024Jun04_A3SS_voom"
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"SE_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by TAM elements, followed by AML elements
piece_exp <- piece[,c(1,3,5,7,2,4,6,8)]
colnames(piece_exp) <- colnames(piece[,c(1,3,5,7,2,4,6,8)])
rownames(piece_exp) <- rownames(piece)
string_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
piece_exp_filename <- paste(paste(paste(paste(date,"expression_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
write.csv(piece_exp$counts,piece_exp_filename,quote=FALSE)
write.csv(rownames(piece),piece_filename,quote=FALSE,row.names=FALSE)
write.csv(fd[,2],string_filename,quote=FALSE,row.names=FALSE)
violin_plot_filename = piece_exp_filename

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 7 groups
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=7)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]
cluster_5_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 5,]
cluster_6_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 6,]
cluster_7_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 7,]

cluster_1_filename <- paste(paste(date, "cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "cluster_4", sep="_"),"csv",sep=".")
cluster_5_filename <- paste(paste(date, "cluster_5", sep="_"),"csv",sep=".")
cluster_6_filename <- paste(paste(date, "cluster_6", sep="_"),"csv",sep=".")
cluster_7_filename <- paste(paste(date, "cluster_7", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_5_genes$geneSymbol,cluster_5_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_6_genes$geneSymbol,cluster_6_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_7_genes$geneSymbol,cluster_7_filename,quote=FALSE,row.names=FALSE)

# -

vwts <- voomWithQualityWeights(dge, design=design, normalize.method="quantile", plot=TRUE)

vwtsfit <- lmFit(vwts, design, weights = vwts$weights )
# no other weighting at this time.
 #* c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0))

vwtsfit <- eBayes(vwtsfit, trend=TRUE)
de_results <- topTable(vwtsfit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(featureData[lookup,2])
head(featureData[lookup,2])

# +
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 9
adjusted_pvalue_threshold <- 0.05

# Select genes that meet both fold change and adjusted p-value criteria
significant_genes <- de_results[
  abs(de_results$logFC) > fold_change_threshold &
  de_results$adj.P.Val < adjusted_pvalue_threshold,
]
dim(significant_genes)
# -

lookup <- rownames(significant_genes)
significant_expression <- dge[lookup,]

significant_out <- pheatmap(significant_expression, 
                            cluster_rows5=TRUE, 
                            show_rownames=FALSE,
                            cluster_cols=TRUE, 
                            annotation_col=df, 
                            scale="row",
                            clustering_method = "ward.D2",
                            clustering_distance_cols = "minkowski", 
                            clustering_distance_rows = "minkowski" )

featureData[head(rownames(significant_expression),5),2]

top_gene_list <- as.matrix(featureData[rownames(significant_expression),2])
length(top_gene_list)

top_significant_genes <- dge[rownames(significant_genes),]

# +
start=1
stop=length(top_gene_list)
date="2024Jun04_A3SS_vwts"
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by TAM elements, followed by AML elements
piece_exp <- piece[,c(1,3,5,7,2,4,6,8)]
colnames(piece_exp) <- colnames(piece[,c(1,3,5,7,2,4,6,8)])
rownames(piece_exp) <- rownames(piece)
string_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
piece_exp_filename <- paste(paste(paste(paste(date,"expression_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
write.csv(piece_exp$counts,piece_exp_filename,quote=FALSE)
write.csv(rownames(piece),piece_filename,quote=FALSE,row.names=FALSE)
write.csv(fd[,2],string_filename,quote=FALSE,row.names=FALSE)
violin_plot_filename = piece_exp_filename
# -

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
save_pheatmap_pdf(outpiece, "2024Jun04_A3SS_vwts_10fold_top_significant_genes.pdf")

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 7 groups
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=7)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]
cluster_5_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 5,]
cluster_6_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 6,]
cluster_7_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 7,]

cluster_1_filename <- paste(paste(date, "cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "cluster_4", sep="_"),"csv",sep=".")
cluster_5_filename <- paste(paste(date, "cluster_5", sep="_"),"csv",sep=".")
cluster_6_filename <- paste(paste(date, "cluster_6", sep="_"),"csv",sep=".")
cluster_7_filename <- paste(paste(date, "cluster_7", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_5_genes$geneSymbol,cluster_5_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_6_genes$geneSymbol,cluster_6_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_7_genes$geneSymbol,cluster_7_filename,quote=FALSE,row.names=FALSE)

# -

# Make violin plots for each of the IDs which are for each row encoded as GeneSymbol.UniqueJunctionIdentifier

# +
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# Read the CSV file (assuming file path is already defined as violin_plot_filename)
data <- read.csv(violin_plot_filename, stringsAsFactors = FALSE)

# Transform the data from wide to long format
data_long <- melt(data, id.vars = "X", variable.name = "Sample", value.name = "Expression")

# Ensure the Sample column is treated as a character vector
data_long$Sample <- as.character(data_long$Sample)

# Extract individual and state information from Sample column
data_long <- data_long %>%
  mutate(Individual = sapply(strsplit(Sample, "\\."), `[`, 1),
         State = sapply(strsplit(Sample, "\\."), `[`, 2))

# Map state codes to state names
state_mapping <- c("03A" = "TAM", "40A" = "AML")
data_long$State <- state_mapping[data_long$State]

# Calculate the mean differences between AML and TAM for sorting
mean_diffs <- data_long %>%
  group_by(X, State) %>%
  summarise(Mean = mean(Expression), .groups = 'drop') %>%
  pivot_wider(names_from = State, values_from = Mean) %>%
  mutate(Diff = AML - TAM) %>%
  arrange(Diff)

# Store the plots in a list
plots <- list()

# Plot violin plots for each gene symbol and splicing junction identifier in sorted order
for (gene_id in mean_diffs$X) {
  gene_data <- subset(data_long, X == gene_id)
  gene_data$State <- factor(gene_data$State, levels = c("TAM", "AML")) # Ensure TAM is plotted first
  
  mean_TAM <- mean(gene_data$Expression[gene_data$State == "TAM"])
  mean_AML <- mean(gene_data$Expression[gene_data$State == "AML"])
  line_data <- data.frame(State = c("TAM", "AML"), Mean = c(mean_TAM, mean_AML))
  
  p <- ggplot(gene_data, aes(x = State, y = Expression, fill = State)) +
    geom_violin() +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) + # Add individual expression values as dots
    geom_line(aes(group = Individual), color = "blue", alpha = 0.5) + # Line connecting samples of the same individual
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) + # Plot mean as points
    ggtitle(paste("Violin Plot for", gene_id)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("TAM" = "skyblue", "AML" = "salmon")) + # Custom colors for states
    geom_line(data = line_data, aes(x = State, y = Mean, group = 1), color = "red", linewidth = 1) # Line connecting means
  
  print(p)  # Print the plot to the notebook
  
  # Save the plot to the list
  plots[[gene_id]] <- p
}
# Load necessary libraries
library(pdftools)

# Directory to save individual PDF files
dir.create("violin_plots")

# Save each plot to a separate PDF file
pdf_files <- c()
for (i in seq_along(plots)) {
  plot_name <- paste0("violin_plots/A3SS_violin_plot_", names(plots)[i], ".pdf")
  ggsave(plot_name, plot = plots[[i]], width = 10, height = 8)
  pdf_files <- c(pdf_files, plot_name)
}

# Combine individual PDF files into a single PDF
pdf_combine(input = pdf_files, output = "A3SS_violin_plots_sorted.pdf")

# Cleanup: remove the individual PDF files
#file.remove(pdf_files)
# -


