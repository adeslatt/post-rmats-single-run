# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
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


BiocManager::install("dplyr")

library(Glimma)
library(dplyr)
library(edgeR)

setwd("/Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/projects/post-rmats-single-run/TAM.AML.all/A3SS_calculate")

getwd()


cts <- as.matrix(read.csv("A3SS.IJC.w.coordinates.matrix.csv",sep=",",row.names="ID"))

dim(cts)

cts[1:3,11:dim(cts)[2]]

featureData <- data.frame(cts[,1:10])
featureData[1:3,]

featureData <- featureData[,c(1,2)]

head(featureData,2)

cts <- data.matrix(cts[,11:dim(cts)[2]])
mode(cts) <- "integer"
is.integer(cts)

dim(cts)
head(cts,2)

colnames(cts)

colnames(cts[,c(113:116)])

# I tried to rename colnames but it did not rename the existing object, rather it created a new object.
#
# From code co-pilot: The issue with your code arises from the fact that when you use subsetting in R with cts[, c(113:116)], it returns a new object, and changing the column names of this new object does not affect the original data frame. To reassign the column names in the original data frame, you need to directly access and modify the colnames attribute of the cts data frame.
#

# +
# Assign new column names directly to the specified columns
colnames(cts)[c(113:116)] <- c("PAUTLA.03A", "PAUTLA.40A", "PAVUDU.03A", "PAVUDU.40A")

# Verify the column names
colnames(cts)[c(113:116)]
# -

# The PAWHSD samples are not TAM and AML but in fact TAM and TAM - the resulting heatmap when included showed they clustered together.  We will eliminate them from subsequent analyses.   To find and replace we need to work on the existing object. Get the mask object (which columns equal PAWHSD).

# +
dim (cts)
# Find column numbers for column names that contain "PAWHSD"
column_numbers <- grep("PAWHSD", colnames(cts))

# Remove these columns from the cts data frame
cts <- cts[, -column_numbers]
dim(cts)
# -

coldata <- read.csv("/Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/projects/post-rmats-single-run/design/design_all_matrix.csv", row.names=1)

head(coldata)

# Create a new column 'name' with rownames
coldata$name <- rownames(coldata)
head(coldata)

# +
# Filter for preAML paired samples
preAML_paired_samples <- coldata %>%
  filter(condition == "preAML" & paired == "yes") %>%
  pull(name)

# Filter for AML paired samples
AML_paired_samples <- coldata %>%
  filter(condition == "AML" & paired == "yes") %>%
  pull(name)

# Filter the unpaired samples
noAML_unpaired_samples <- coldata %>% 
   filter(condition == "noAML" & paired == "no") %>%
   pull (name)

preAML_unpaired_samples <- coldata %>%
   filter (condition == "preAML" & paired == "no") %>%
   pull (name)
length(preAML_paired_samples)
length(AML_paired_samples)
length(preAML_unpaired_samples)
length(noAML_unpaired_samples)

# +
# Select columns corresponding to preAML paired samples
cts_preAML_paired <- cts[, colnames(cts) %in% preAML_paired_samples]

# Select columns corresponding to AML paired samples
cts_AML_paired <- cts[, colnames(cts) %in% AML_paired_samples]

# Select columns corresponding to TAM (no progress to AML) unpaired samples
cts_noAML_unpaired <- cts[, colnames(cts) %in% noAML_unpaired_samples]

# Select columns corresponding to TAM (no progress to AML) unpaired samples
cts_preAML_unpaired <- cts[, colnames(cts) %in% preAML_unpaired_samples]

# Display the first few rows of the filtered dataframes
head(cts_preAML_paired)
head(cts_AML_paired)
head(cts_noAML_unpaired)
head(cts_preAML_unpaired)


# +

# Ensure colnames of cts match the 'name' column in design_all_matrix
# Assuming 'cts' is already loaded in the R environment
#matching_cols <- colnames(cts) %in% 

all(coldata$name %in% colnames(cts))


# +
missing_names <- coldata$name[!coldata$name %in% colnames(cts)]

# Print the missing names
cat("The following names are missing from cts:\n")
print(missing_names)

# -

# This is an interesting catch - this was an inadvertant add -- I had removed it from all the references so I don't include it going forward.  But did not remove it from the design matrix.

# **Key Points:**
#
# 1. **Threshold Definition:** Define a count threshold (e.g., 1000).
# 2. **Criteria for Selection:** Define a function to select rows that have counts above the threshold in at least half of the samples in each condition.
# 3. **Row Selection:** Apply this criteria to each condition's sub-matrix.
# 4. **Combine Rows:** Combine the selected rows from all conditions to form the final matrix.
# 5. **MDS Plot:** Proceed with the MDS plot using the final matrix.
# This approach ensures that the rows included in the final matrix have significant counts in most of the samples for each condition, which should help to maximize the differences in the signals.

# Select columns corresponding to preAML paired samples
cts_paired <- cts[, c(colnames(cts) %in% preAML_paired_samples | colnames(cts) %in% AML_paired_samples)]
dim(cts_paired)


cts_unpaired <- cts[, c(colnames(cts) %in% preAML_unpaired_samples | colnames(cts) %in% noAML_unpaired_samples)]
dim(cts_unpaired)


colnames(cts_paired)

coldata_paired <- read.csv("/Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/projects/post-rmats-single-run/design/design_paired_matrix.csv", row.names=1)
coldata_unpaired_noAML_paired_AML  <- read.csv("/Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/projects/post-rmats-single-run/design/design_unpaired_noAML_paired_AML.csv", row.names=1)
coldata_unpaired_preAML_paired_AML <- read.csv("/Users/annedeslattesmays/Scitechcon\ Dropbox/Anne\ DeslattesMays/projects/post-rmats-single-run/design/design_unpaired_preAML_paired_AML.csv", row.names=1)


head(coldata)
head(coldata_unpaired_noAML_paired_AML)
head(coldata_unpaired_preAML_paired_AML)

coldata_paired
dim(coldata_paired)

# +
# Load necessary libraries
library(dplyr)

# Assuming coldata is already defined and has row.names
# Assuming cts_paired is already defined

# Identify the columns for preAML and paired == "yes"
paired_preAML_cols <- rownames(coldata)[coldata$condition == "preAML" & coldata$paired == "yes"]

# Identify the columns for AML and paired == "yes"
paired_AML_cols <- rownames(coldata)[coldata$condition == "AML" & coldata$paired == "yes"]

# Convert column names to indices
paired_preAML_indices <- which(colnames(cts_paired) %in% paired_preAML_cols)
paired_AML_indices <- which(colnames(cts_paired) %in% paired_AML_cols)

# Define the threshold
threshold <- 10  # You can change this to 10, 100, etc.

# Condition 1: Rows with count > threshold in all paired preAML columns
paired_preAML_rows_condition <- apply(cts_paired[, paired_preAML_indices], 1, function(row) all(row > threshold))
paired_preAML_matrix <- cts_paired[paired_preAML_rows_condition, paired_preAML_indices]
dim(paired_preAML_matrix)

# Condition 2: Rows with count > threshold in all paired AML columns
paired_AML_rows_condition <- apply(cts_paired[, paired_AML_indices], 1, function(row) all(row > threshold))
paired_AML_matrix <- cts_paired[paired_AML_rows_condition, paired_AML_indices]
dim(paired_AML_matrix)


# +

# Identify the columns for preAML and paired == "no"
unpaired_preAML_cols <- rownames(coldata)[coldata$condition == "preAML" & coldata$paired == "no"]

# Identify the columns for noAML and paired == "no"
unpaired_noAML_cols <- rownames(coldata)[coldata$condition == "noAML" & coldata$paired == "no"]

# Convert column names to indices
unpaired_preAML_indices <- which(colnames(cts_unpaired) %in% unpaired_preAML_cols)
unpaired_noAML_indices  <- which(colnames(cts_unpaired) %in% unpaired_noAML_cols)

# Condition 1: Rows with count > threshold in all unpaired preAML columns
unpaired_preAML_rows_condition <- apply(cts_unpaired[, unpaired_preAML_indices], 1, function(row) all(row > threshold))
unpaired_preAML_matrix <- cts_unpaired[unpaired_preAML_rows_condition, unpaired_preAML_indices]
dim(unpaired_preAML_matrix)

# Condition 2: Rows with count > threshold in all unpaired noAML columns
unpaired_noAML_rows_condition <- apply(cts_unpaired[, unpaired_noAML_indices], 1, function(row) all(row > threshold))
unpaired_noAML_matrix <- cts_unpaired[unpaired_noAML_rows_condition, unpaired_noAML_indices]
dim(unpaired_noAML_matrix)


# +
# Extract sample names for different conditions
paired_preAML_samples <- coldata$name[coldata$condition == "preAML" & coldata$paired == "yes"]
paired_AML_samples    <- coldata$name[coldata$condition == "AML"    & coldata$paired == "yes"]

print(paired_preAML_samples)
print(paired_AML_samples)

# Filter paired data based on the threshold
paired_preAML_columns <- colnames(cts) %in% paired_preAML_samples
paired_AML_columns    <- colnames(cts) %in% paired_AML_samples

# Select columns for paired and unpaired samples
paired_data <- cts[, colnames(cts) %in% c(paired_preAML_samples, paired_AML_samples)]
dim(paired_data)
colnames(paired_data)
#
# make the final paired_matrix
#
# Combine the sub-matrices by keeping rows that satisfy either condition 1 or condition 2
#
paired_final_matrix <- cts_paired[paired_AML_rows_condition | paired_preAML_rows_condition, ]
dim(paired_final_matrix)
head(paired_final_matrix)

# -
#
# now make one of two unpaired preAML matrices
#
# 1. unpaired_preAML and paired AML
# 2. unpaired_noAML and paired AML
#

unpaired_preAML_samples <- coldata$name[coldata$condition == "preAML" & coldata$paired == "no"]
unpaired_noAML_samples  <- coldata$name[coldata$condition == "noAML" & coldata$paired == "no"]
unpaired_final_matrix   <- cts_unpaired[unpaired_noAML_rows_condition | unpaired_preAML_rows_condition, ]
dim(unpaired_final_matrix)

# +
# can we make two matrices 
# 1. combining unpaired preAML with paired AML
# and a matrix
# 2. combining unpaired noAML with paired AML
#
combined_unpaired_preAML_paired_AML_data <- cts[, colnames(cts) %in% c(unpaired_preAML_samples, paired_AML_samples)]
combined_unpaired_noAML_paired_AML_data  <- cts[, colnames(cts) %in% c(unpaired_noAML_samples, paired_AML_samples)]
combined_final_unpaired_preAML_paired_AML_matrix <- combined_unpaired_preAML_paired_AML_data[unpaired_preAML_rows_condition | 
                                                                                             paired_AML_rows_condition,]
combined_final_unpaired_noAML_paired_AML_matrix  <- combined_unpaired_noAML_paired_AML_data [unpaired_noAML_rows_condition | 
                                                                                             paired_AML_rows_condition,]

dim(combined_final_unpaired_preAML_paired_AML_matrix)
dim(combined_final_unpaired_noAML_paired_AML_matrix)


coldata_unpaired_preAML_paired_AML$condition <- factor (coldata_unpaired_preAML_paired_AML$condition)
coldata_unpaired_noAML_paired_AML$condition <- factor (coldata_unpaired_noAML_paired_AML$condition)

levels(coldata_unpaired_noAML_paired_AML$condition)
levels(coldata_unpaired_preAML_paired_AML$condition)

# +


# +
# Load necessary libraries
library(edgeR)
library(limma)


#
# Three limma/voom analyses to run:
# 1. paired analyses preAML to AML
# 2. unpaired analyses unpaired preAML to paired AML
# 3. unpaired analyses unpaired noAML to paired AML
#

# 1. Paired analyses
#    key input matrices are:
#    a. coldata_paired
#    b. paired_final_matrix
#
# Assuming paired_final_matrix and coldata_paired are already defined
paired_cts <- paired_final_matrix

# Check dimensions
print(dim(paired_cts))
print(dim(coldata_paired))

# Ensure coldata_paired has the same number of rows as paired_cts has columns
if (nrow(coldata_paired) != ncol(paired_cts)) {
  stop("The number of samples in coldata_paired and paired_cts do not match.")
}

#
# Begin Differential count analyses
# 
# Convert the count data to a DGEList object
#
dge <- DGEList(counts = paired_cts)

# Normalize the data using the TMM method
dge <- calcNormFactors(dge)

# Compute the log-transformed counts per million (CPM)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# Ensure coldata_paired$condition is a factor and create a color palette
coldata_paired$condition <- factor(coldata_paired$condition, levels = c("preAML", "AML"))
color_palette <- c("blue", "green")

# Map the conditions to their respective colors
condition_colors <- color_palette[as.numeric(coldata_paired$condition)]

# Check lengths of logCPM and condition_colors
print(dim(logCPM))
print(length(condition_colors))

# Ensure the lengths match
if (ncol(logCPM) != length(condition_colors)) {
  stop("The number of col in logCPM and the length of condition_colors do not match.")
}

# Create MDS plot with colors for each condition
plotMDS(logCPM, col = condition_colors, labels = coldata_paired$condition)
legend("topright", legend = levels(coldata_paired$condition), col = color_palette, pch = 16)

group <- factor(coldata_paired$condition)

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

de_results <- topTable(fit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(lookup)
head(featureData[lookup,2])

# +
# There are too many values - lets reduce the size a bit more
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 9
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

# voom 
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
date="2024July30_A3SS_paired_preAML_paired_AML_voom"
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by preAML elements, followed by AML elements
piece_exp <- piece[,c((colnames(cts) %in% preAML_paired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% preAML_paired_samples),(colnames(cts) %in% AML_paired_samples))])
rownames(piece_exp) <- rownames(piece)
string_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
piece_exp_filename <- paste(paste(paste(paste(date,"_expression_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
write.csv(piece_exp$counts,piece_exp_filename,quote=FALSE)
write.csv(rownames(piece),piece_filename,quote=FALSE,row.names=FALSE)
write.csv(fd[,2],string_filename,quote=FALSE,row.names=FALSE)

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 major groups
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)

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
head(significant_genes)


head(featureData[lookup,2])


head(significant_expression)


# +
paired_preAML_paired_AML_final_matrix <- significant_expression$counts

rownames(paired_preAML_paired_AML_final_matrix) <- lookup
head(paired_preAML_paired_AML_final_matrix)
write.csv(paired_preAML_paired_AML_final_matrix, "paired_preAML_paired_AML_final_matrix.csv",quote=FALSE)
# -

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
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by preAML elements, followed by AML elements
piece_exp <- piece[,c((colnames(cts) %in% preAML_paired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% preAML_paired_samples),(colnames(cts) %in% AML_paired_samples))])
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
save_pheatmap_pdf(outpiece, paste(date,"final_matrix_vwts_10fold_top_significant_genes.pdf",sep="_")

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 groups - should check the dendogram every time to choose major grouping
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_vwts_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_vwts_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_vwts_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_vwts_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)

#
# 2. unpaired analyses unpaired preAML to paired AML
#    key input matrices are:
#    a. coldata_unpaired_preAML_paired_AML
#    b. combined_final_unpaired_preAML_paired_AML_cts
#

# +

combined_final_unpaired_preAML_paired_AML_cts <- combined_final_unpaired_preAML_paired_AML_matrix

# Check dimensions
print(dim(combined_final_unpaired_preAML_paired_AML_cts))
print(dim(coldata_unpaired_preAML_paired_AML))

# Ensure matrix has the same number of rows as cts has columns
if (nrow(coldata_unpaired_preAML_paired_AML) != ncol(combined_final_unpaired_preAML_paired_AML_cts)) {
  stop("The number of samples in coldata_paired and paired_cts do not match.")
}

# Convert the count data to a DGEList object
dge <- DGEList(counts = combined_final_unpaired_preAML_paired_AML_cts)

# Normalize the data using the TMM method
dge <- calcNormFactors(dge)

# Compute the log-transformed counts per million (CPM)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# Ensure coldata_paired$condition is a factor and create a color palette
coldata_unpaired_preAML_paired_AML$condition <- factor(coldata_unpaired_preAML_paired_AML$condition, levels = c("preAML", "AML"))
color_palette <- c("blue", "green")

# Map the conditions to their respective colors
condition_colors <- color_palette[as.numeric(coldata_unpaired_preAML_paired_AML$condition)]

# Check lengths of logCPM and condition_colors
print(dim(logCPM))
print(length(condition_colors))

# Ensure the lengths match
if (ncol(logCPM) != length(condition_colors)) {
  stop("The number of col in logCPM and the length of condition_colors do not match.")
}

# Create MDS plot with colors for each condition
plotMDS(logCPM, col = condition_colors, labels = coldata_unpaired_preAML_paired_AML$condition)
legend("topright", legend = levels(coldata_unpaired_preAML_paired_AML$condition), col = color_palette, pch = 16)


# +
coldata_unpaired_preAML_paired_AML
# Create a factor for the conditions
group <- factor(coldata_unpaired_preAML_paired_AML$condition)

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

de_results <- topTable(fit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(lookup)
head(featureData[lookup,2])

# +
# There are too many values - lets reduce the size a bit more
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 9
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

# voom 
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
date="2024July30_A3SS_unpaired_preAML_paired_AML_voom"
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by preAML elements, followed by AML element
piece_exp <- piece[,c((colnames(cts) %in% preAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% preAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))])

rownames(piece_exp) <- rownames(piece)
string_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
piece_exp_filename <- paste(paste(paste(paste(date,"_expression_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
write.csv(piece_exp$counts,piece_exp_filename,quote=FALSE)
write.csv(rownames(piece),piece_filename,quote=FALSE,row.names=FALSE)
write.csv(fd[,2],string_filename,quote=FALSE,row.names=FALSE)

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 major groups
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)

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
head(significant_genes)


head(featureData[lookup,2])


head(significant_expression)


# +
unpaired_preAML_paired_AML_final_matrix <- significant_expression$counts

rownames(unpaired_preAML_paired_AML_final_matrix) <- lookup
head(unpaired_preAML_paired_AML_final_matrix)
write.csv(unpaired_preAML_paired_AML_final_matrix, "unpaired_preAML_paired_AML_final_matrix.csv",quote=FALSE)
# -

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
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by TAM elements, followed by AML elements
piece_exp <- piece[,c((colnames(cts) %in% preAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% preAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))])
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
save_pheatmap_pdf(outpiece, "2024Jul30_A3SS_unpaired_preAML_paired_AML_final_matrix_vwts_10fold_top_significant_genes.pdf")

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 groups - should check the dendogram every time to choose major grouping
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)

#
# 3. unpaired analyses unpaired noAML to paired AML
#    key input matrices are:
#    a. coldata_unpaired_noAML_paired_AML
#    b. combined_final_unpaired_noAML_paired_AML_cts
#

# +

combined_final_unpaired_noAML_paired_AML_cts <- combined_final_unpaired_noAML_paired_AML_matrix

# Check dimensions
print(dim(combined_final_unpaired_noAML_paired_AML_cts))
print(dim(coldata_unpaired_noAML_paired_AML))

# Ensure matrix has the same number of rows as cts has columns
if (nrow(coldata_unpaired_noAML_paired_AML) != ncol(combined_final_unpaired_noAML_paired_AML_cts)) {
  stop("The number of samples in coldata_paired and paired_cts do not match.")
}

# Convert the count data to a DGEList object
dge <- DGEList(counts = combined_final_unpaired_noAML_paired_AML_cts)

# Normalize the data using the TMM method
dge <- calcNormFactors(dge)

# Compute the log-transformed counts per million (CPM)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# Ensure coldata_paired$condition is a factor and create a color palette
coldata_unpaired_noAML_paired_AML$condition <- factor(coldata_unpaired_noAML_paired_AML$condition, levels = c("noAML", "AML"))
color_palette <- c("blue", "green")

# Map the conditions to their respective colors
condition_colors <- color_palette[as.numeric(coldata_unpaired_noAML_paired_AML$condition)]

# Check lengths of logCPM and condition_colors
print(dim(logCPM))
print(length(condition_colors))

# Ensure the lengths match
if (ncol(logCPM) != length(condition_colors)) {
  stop("The number of col in logCPM and the length of condition_colors do not match.")
}

# Create MDS plot with colors for each condition
plotMDS(logCPM, col = condition_colors, labels = coldata_unpaired_noAML_paired_AML$condition)
legend("topright", legend = levels(coldata_unpaired_noAML_paired_AML$condition), col = color_palette, pch = 16)


# +
coldata_unpaired_noAML_paired_AML
# Create a factor for the conditions
group <- factor(coldata_unpaired_noAML_paired_AML$condition)

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

de_results <- topTable(fit, coef=ncol(design), n=Inf) 
lookup <- rownames(de_results)
length(lookup)
head(featureData[lookup,2])

# +
# There are too many values - lets reduce the size a bit more
# Assuming you have the 'de_results' object from topTable
fold_change_threshold <- 9
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

# voom 
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
date="2024July30_A3SS_unpaired_noAML_paired_AML_voom"
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by noAML elements, followed by AML element
piece_exp <- piece[,c((colnames(cts) %in% noAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% noAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))])

rownames(piece_exp) <- rownames(piece)
string_filename <- paste(paste(paste(paste(date,"_string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
piece_exp_filename <- paste(paste(paste(paste(date,"_expression_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")
write.csv(piece_exp$counts,piece_exp_filename,quote=FALSE)
write.csv(rownames(piece),piece_filename,quote=FALSE,row.names=FALSE)
write.csv(fd[,2],string_filename,quote=FALSE,row.names=FALSE)

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 major groups
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)

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
head(significant_genes)


head(featureData[lookup,2])


head(significant_expression)


# +
unpaired_noAML_paired_AML_final_matrix <- significant_expression$counts

rownames(unpaired_noAML_paired_AML_final_matrix) <- lookup
head(unpaired_noAML_paired_AML_final_matrix)
write.csv(unpaired_noAML_paired_AML_final_matrix, "unpaired_noAML_paired_AML_final_matrix.csv",quote=FALSE)
# -

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
piece <-top_significant_genes[significant_out$tree_row$order[start:stop],]
fd <- data.frame(featureData[rownames(piece),])
genejunction <- paste(featureData[rownames(piece),2],rownames(piece),sep=".")
rownames(fd) <- genejunction
rownames(piece) <- genejunction
outpiece<-pheatmap(piece, cluster_rows5=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_cols = "minkowski", clustering_distance_rows = "minkowski" )
piece_filename <- paste(paste(paste(paste(date,"string_top_gene_list",sep="_"),start,sep="_"),stop,sep="_"),"csv",sep=".")

# Order by TAM elements, followed by AML elements
piece_exp <- piece[,c((colnames(cts) %in% noAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))]
colnames(piece_exp) <- colnames(piece[,c((colnames(cts) %in% noAML_unpaired_samples),(colnames(cts) %in% AML_paired_samples))])
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
save_pheatmap_pdf(outpiece, "2024Jul30_A3SS_unpaired_noAML_paired_AML_final_matrix_vwts_10fold_top_significant_genes.pdf")

# +
#If you want something like gene-to-cluster assignment, you can 'cut' your row dendrogram into a pre-selected number of groups as follows:
# -- inspecating above the rows seem to fall into 4 groups - should check the dendogram every time to choose major grouping
clusters<- as.matrix(row_clusters<- sort(cutree(significant_out$tree_row, k=4)),nrows=dim(top_genes_expression)[1],ncols=1)
genes_in_clusters = featureData[rownames(clusters),2]
genes_in_clusters.df <- data.frame(featureData[rownames(clusters),2], clusters)
colnames(genes_in_clusters.df) <- c("geneSymbol","cluster")
dim(genes_in_clusters.df)
cluster_1_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 1,]
cluster_2_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 2,]
cluster_3_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 3,]
cluster_4_genes <- genes_in_clusters.df[genes_in_clusters.df$cluster == 4,]

cluster_1_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_1", sep="_"),"csv",sep=".")
cluster_2_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_2", sep="_"),"csv",sep=".")
cluster_3_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_3", sep="_"),"csv",sep=".")
cluster_4_filename <- paste(paste(date, "_paired_A3SS_vwts_cluster_4", sep="_"),"csv",sep=".")

write.csv(cluster_1_genes$geneSymbol,cluster_1_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_2_genes$geneSymbol,cluster_2_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_3_genes$geneSymbol,cluster_3_filename,quote=FALSE,row.names=FALSE)
write.csv(cluster_4_genes$geneSymbol,cluster_4_filename,quote=FALSE,row.names=FALSE)



# -

# Make violin plots for each of the IDs which are for each row encoded as GeneSymbol.UniqueIdentifier

# +
# Load necessary libraries
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(pdftools)
library(gridExtra)

# Read the unpaired and paired datasets
unpaired_data <- read.csv("unpaired_final_matrix.csv", stringsAsFactors = FALSE)
paired_data <- read.csv("paired_final_matrix.csv", stringsAsFactors = FALSE)

# Transform paired data from wide to long format
paired_data_long <- melt(paired_data, variable.name = "Sample", value.name = "Expression", id.vars = "X")
paired_data_long <- paired_data_long %>%
  mutate(Individual = sapply(strsplit(as.character(Sample), "\\."), `[`, 1),
         State = sapply(strsplit(as.character(Sample), "\\."), `[`, 2))

# Transform unpaired data from wide to long format
unpaired_data_long <- melt(unpaired_data, variable.name = "Sample", value.name = "Expression", id.vars = "X")

# Ensure the Sample column is treated as a character vector
unpaired_data_long$Sample <- as.character(unpaired_data_long$Sample)

# Extract the state from the Sample column and create the Individual column
unpaired_data_long <- unpaired_data_long %>%
  mutate(Individual = sapply(strsplit(Sample, "\\."), `[`, 1))

# Map sample names to conditions using the coldata DataFrame
unpaired_data_long <- unpaired_data_long %>%
  left_join(coldata, by = c("Sample" = "name")) %>%
  mutate(State = condition) %>%
  select(-condition)

# Create directories
dir.create("paired_violin_plots", showWarnings = FALSE)
dir.create("comparison_violin_plots", showWarnings = FALSE)
dir.create("comparison_violin_plots/consistent", showWarnings = FALSE)
dir.create("comparison_violin_plots/opposite", showWarnings = FALSE)

# -

dim(unpaired_data)
dim(paired_data)



# Function to create paired violin plots
create_paired_violin_plot <- function(gene_data, gene_name, filename, t_test_mean, t_test_median) {
  gene_data$State <- factor(gene_data$State, levels = c("preAML", "AML"))
  
  mean_preAML <- mean(gene_data$Expression[gene_data$State == "preAML"])
  mean_AML <- mean(gene_data$Expression[gene_data$State == "AML"])
  median_preAML <- median(gene_data$Expression[gene_data$State == "preAML"])
  median_AML <- median(gene_data$Expression[gene_data$State == "AML"])
  
  line_data <- data.frame(State = c("preAML", "AML"), Mean = c(mean_preAML, mean_AML), Median = c(median_preAML, median_AML))
  
  p <- ggplot(gene_data, aes(x = State, y = Expression, fill = State)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) + # Add individual expression values as dots
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) + # Plot mean as points
    geom_line(data = line_data, aes(x = State, y = Mean, group = 1), color = "blue", size = 1) + # Line connecting means
    geom_line(data = line_data, aes(x = State, y = Median, group = 1), color = "green", size = 1) + # Line connecting medians
    geom_line(aes(group = Individual), color = "grey", alpha = 0.5) + # Line connecting paired samples
    ggtitle(paste("Paired Violin Plot for", gene_name, "\nMean p-value: ", signif(t_test_mean$p.value, 3), " Median p-value: ", signif(t_test_median$p.value, 3))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("preAML" = "skyblue", "AML" = "salmon")) # Custom colors for states
  
  ggsave(filename, plot = p, width = 10, height = 8)
  return(p)
}


# Function to create comparison violin plots
create_comparison_violin_plot <- function(gene_data, gene_name, filename, left_state, right_state, t_test_mean, t_test_median) {
  gene_data$State <- factor(gene_data$State, levels = c(left_state, right_state))
  
  mean_left <- mean(gene_data$Expression[gene_data$State == left_state])
  mean_right <- mean(gene_data$Expression[gene_data$State == right_state])
  median_left <- median(gene_data$Expression[gene_data$State == left_state])
  median_right <- median(gene_data$Expression[gene_data$State == right_state])
  
  line_data <- data.frame(State = c(left_state, right_state), Mean = c(mean_left, mean_right), Median = c(median_left, median_right))
  
  p <- ggplot(gene_data, aes(x = State, y = Expression, fill = State)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) + # Add individual expression values as dots
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) + # Plot mean as points
    geom_line(data = line_data, aes(x = State, y = Mean, group = 1), color = "blue", size = 1) + # Line connecting means
    geom_line(data = line_data, aes(x = State, y = Median, group = 1), color = "green", size = 1) + # Line connecting medians
    ggtitle(paste("Violin Plot for", gene_name, "\nMean p-value: ", signif(t_test_mean$p.value, 3), " Median p-value: ", signif(t_test_median$p.value, 3))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c(left_state = "skyblue", right_state = "salmon")) # Custom colors for states
  
  ggsave(filename, plot = p, width = 10, height = 8)
  return(p)
}

# +
# Initialize lists to store statistical results
stat_results_consistent <- list()
stat_results_opposite <- list()

# Plot paired data
for (gene_id in unique(paired_data_long$X)) {
  gene_name <- featureData[gene_id, 2]
  gene_data <- paired_data_long %>% filter(X == gene_id)
  
  # Perform t-tests between means
  t_test_mean <- t.test(gene_data$Expression[gene_data$State == "preAML"], 
                        gene_data$Expression[gene_data$State == "AML"])
  
  # Perform Wilcoxon tests between medians
  t_test_median <- wilcox.test(gene_data$Expression[gene_data$State == "preAML"], 
                               gene_data$Expression[gene_data$State == "AML"])
  
  create_paired_violin_plot(gene_data, gene_name, paste0("paired_violin_plots/", gene_name, "_paired_violin_plot.pdf"), t_test_mean, t_test_median)
}


# +
# Plot consistent test: Paired AML with unpaired preAML
consistent_plots <- list()
for (gene_id in unique(paired_data_long$X)) {
  gene_name <- featureData[gene_id, 2]
  
  paired_gene_data <- paired_data_long %>% filter(X == gene_id & State == "AML") %>% mutate(Source = "Paired")
  unpaired_gene_data <- unpaired_data_long %>% filter(X == gene_id & State == "preAML") %>% mutate(State = "preAML_unpaired", Source = "Unpaired")
  
  if (nrow(unpaired_gene_data) > 0) {
    # Ensure columns match before combining
    common_cols <- intersect(names(paired_gene_data), names(unpaired_gene_data))
    gene_data <- rbind(paired_gene_data[, common_cols], unpaired_gene_data[, common_cols])
  
    # Perform t-tests between means
    t_test_mean <- t.test(gene_data$Expression[gene_data$State == "preAML_unpaired"], 
                          gene_data$Expression[gene_data$State == "AML"])
    
    # Perform Wilcoxon tests between medians
    t_test_median <- wilcox.test(gene_data$Expression[gene_data$State == "preAML_unpaired"], 
                                 gene_data$Expression[gene_data$State == "AML"])
    
    p <- create_comparison_violin_plot(gene_data, gene_name, paste0("comparison_violin_plots/consistent/", gene_name, "_consistent_violin_plot.pdf"), "preAML_unpaired", "AML", t_test_mean, t_test_median)
    consistent_plots[[gene_id]] <- p
    
    # Store statistical results
    stat_results_consistent[[gene_id]] <- list(
      gene_name = gene_name,
      mean_p_value = t_test_mean$p.value,
      median_p_value = t_test_median$p.value
    )
  }
}

# Plot opposite test: Paired AML with unpaired noAML
opposite_plots <- list()
for (gene_id in unique(paired_data_long$X)) {
  gene_name <- featureData[gene_id, 2]
  
  paired_gene_data <- paired_data_long %>% filter(X == gene_id & State == "AML") %>% mutate(Source = "Paired")
  unpaired_gene_data <- unpaired_data_long %>% filter(X == gene_id & State == "noAML") %>% mutate(State = "noAML", Source = "Unpaired")
  
  if (nrow(unpaired_gene_data) > 0) {
    # Ensure columns match before combining
    common_cols <- intersect(names(paired_gene_data), names(unpaired_gene_data))
    gene_data <- rbind(paired_gene_data[, common_cols], unpaired_gene_data[, common_cols])
  
    # Perform t-tests between means
    t_test_mean <- t.test(gene_data$Expression[gene_data$State == "noAML"], 
                          gene_data$Expression[gene_data$State == "AML"])
    
    # Perform Wilcoxon tests between medians
    t_test_median <- wilcox.test(gene_data$Expression[gene_data$State == "noAML"], 
                                 gene_data$Expression[gene_data$State == "AML"])
    
    p <- create_comparison_violin_plot(gene_data, gene_name, paste0("comparison_violin_plots/opposite/", gene_name, "_opposite_violin_plot.pdf"), "noAML", "AML", t_test_mean, t_test_median)
    opposite_plots[[gene_id]] <- p
    
    # Store statistical results
    stat_results_opposite[[gene_id]] <- list(
      gene_name = gene_name,
      mean_p_value = t_test_mean$p.value,
      median_p_value = t_test_median$p.value
    )
  }
}

# Combine statistical results into data frames
stat_results_consistent_df <- do.call(rbind, lapply(stat_results_consistent, function(x) data.frame(gene_name = x$gene_name, mean_p_value = x$mean_p_value, median_p_value = x$median_p_value)))
stat_results_opposite_df <- do.call(rbind, lapply(stat_results_opposite, function(x) data.frame(gene_name = x$gene_name, mean_p_value = x$mean_p_value, median_p_value = x$median_p_value)))

# Save statistical results to CSV files
write.csv(stat_results_consistent_df, "stat_results_consistent.csv", row.names = FALSE)
write.csv(stat_results_opposite_df, "stat_results_opposite.csv", row.names = FALSE)

# Display the plots (optional, can be removed if running in a non-interactive environment)
for (p in consistent_plots) {
  print(p)
}
for (p in opposite_plots) {
  print(p)
}

# Combine individual PDF files into a single PDF for consistent and opposite plots
pdf_combine(input = list.files("comparison_violin_plots/consistent", full.names = TRUE), output = "consistent_violin_plots_combined.pdf")
pdf_combine(input = list.files("comparison_violin_plots/opposite", full.names = TRUE), output = "opposite_violin_plots_combined.pdf")

# +
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Calculate summary statistics for opposite test p-values
opposite_p_values <- stat_results_opposite_df %>%
  mutate(Significant = ifelse(mean_p_value < 0.05, "Yes", "No"))

summary_stats_opposite <- opposite_p_values %>%
  summarize(
    mean_p_value_mean = mean(mean_p_value),
    median_p_value_mean = median(mean_p_value),
    sd_p_value_mean = sd(mean_p_value),
    mean_p_value_median = mean(median_p_value),
    median_p_value_median = median(median_p_value),
    sd_p_value_median = sd(median_p_value)
  )

print(summary_stats_opposite)

# Visualize the distribution of mean p-values
ggplot(opposite_p_values, aes(x = mean_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Mean P-values for Opposite Test", x = "Mean P-value", y = "Frequency")

# Visualize the distribution of median p-values
ggplot(opposite_p_values, aes(x = median_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Median P-values for Opposite Test", x = "Median P-value", y = "Frequency")

# QQ plot for mean p-values
ggplot(opposite_p_values, aes(sample = mean_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Mean P-values for Opposite Test")

# QQ plot for median p-values
ggplot(opposite_p_values, aes(sample = median_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Median P-values for Opposite Test")

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(gene_name, mean_p_value)

# Identify significant genes based on median p-values
significant_genes_median_opposite <- opposite_p_values %>%
  filter(median_p_value < 0.05) %>%
  select(gene_name, median_p_value)

# Print the number of significant genes
cat("Number of significant genes based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")
cat("Number of significant genes based on median p-value (opposite test):", nrow(significant_genes_median_opposite), "\n")


# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)
write.csv(significant_genes_median_opposite, "significant_genes_median_opposite.csv", row.names = FALSE)
# -



# +
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Calculate summary statistics for opposite test p-values
consistent_p_values <- stat_results_consistent_df %>%
  mutate(Significant = ifelse(mean_p_value < 0.05, "Yes", "No"))

summary_stats_consistent <- consistent_p_values %>%
  summarize(
    mean_p_value_mean = mean(mean_p_value),
    median_p_value_mean = median(mean_p_value),
    sd_p_value_mean = sd(mean_p_value),
    mean_p_value_median = mean(median_p_value),
    median_p_value_median = median(median_p_value),
    sd_p_value_median = sd(median_p_value)
  )

print(summary_stats_consistent)

# Visualize the distribution of mean p-values
ggplot(consistent_p_values, aes(x = mean_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Mean P-values for Consistent Test", x = "Mean P-value", y = "Frequency")

# Visualize the distribution of median p-values
ggplot(consistent_p_values, aes(x = median_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Median P-values for Consistent Test", x = "Median P-value", y = "Frequency")

# QQ plot for mean p-values
ggplot(consistent_p_values, aes(sample = mean_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Mean P-values for Opposite Test")

# QQ plot for median p-values
ggplot(consistent_p_values, aes(sample = median_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Median P-values for Opposite Test")

# Identify significant genes based on mean p-values
significant_genes_mean_consistent <- consistent_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(gene_name, mean_p_value)

# Identify significant genes based on median p-values
significant_genes_median_consistent < consistent_p_values %>%
  filter(median_p_value < 0.05) %>%
  select(gene_name, median_p_value)

# Print the number of significant genes
cat("Number of significant genes based on mean p-value (consistent test):", nrow(significant_genes_mean_consistent), "\n")
cat("Number of significant genes based on median p-value (consistent test):", nrow(significant_genes_median_consistent), "\n")


# Save significant genes to CSV files
write.csv(significant_genes_mean_consistent, "significant_genes_mean_consistent.csv", row.names = FALSE)
write.csv(significant_genes_median_consistent, "significant_genes_median_consistent.csv", row.names = FALSE)

# +
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Calculate summary statistics for opposite test p-values
significant_p_values <- stat_results_opposite_df %>%
  mutate(Significant = ifelse(mean_p_value < 0.05, "Yes", "No"))

summary_stats_opposite <- opposite_p_values %>%
  summarize(
    mean_p_value_mean = mean(mean_p_value),
    median_p_value_mean = median(mean_p_value),
    sd_p_value_mean = sd(mean_p_value),
    mean_p_value_median = mean(median_p_value),
    median_p_value_median = median(median_p_value),
    sd_p_value_median = sd(median_p_value)
  )

print(summary_stats_opposite)

# Visualize the distribution of mean p-values
ggplot(opposite_p_values, aes(x = mean_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Mean P-values for Opposite Test", x = "Mean P-value", y = "Frequency")

# Visualize the distribution of median p-values
ggplot(opposite_p_values, aes(x = median_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Median P-values for Opposite Test", x = "Median P-value", y = "Frequency")

# QQ plot for mean p-values
ggplot(opposite_p_values, aes(sample = mean_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Mean P-values for Opposite Test")

# QQ plot for median p-values
ggplot(opposite_p_values, aes(sample = median_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Median P-values for Opposite Test")

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(gene_name, mean_p_value)

# Identify significant genes based on median p-values
significant_genes_median_opposite <- opposite_p_values %>%
  filter(median_p_value < 0.05) %>%
  select(gene_name, median_p_value)

# Print the number of significant genes
cat("Number of significant genes based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")
cat("Number of significant genes based on median p-value (opposite test):", nrow(significant_genes_median_opposite), "\n")


# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)
write.csv(significant_genes_median_opposite, "significant_genes_median_opposite.csv", row.names = FALSE)

# +
# Combine slopes data for consistent and opposite analysis
slope_data_consistent <- paired_slopes %>%
  left_join(unpaired_slopes %>% select(X, Unpaired_Slope_preAML), by = "X") %>%
  filter(!is.na(Paired_Slope) & !is.na(Unpaired_Slope_preAML)) %>%
  mutate(Consistent = sign(Paired_Slope) == sign(Unpaired_Slope_preAML))

slope_data_opposite <- paired_slopes %>%
  left_join(unpaired_slopes %>% select(X, Unpaired_Slope_noAML), by = "X") %>%
  filter(!is.na(Paired_Slope) & !is.na(Unpaired_Slope_noAML)) %>%
  mutate(Opposite = sign(Paired_Slope) != sign(Unpaired_Slope_noAML))

# Identify significant junctions in consistent and opposite sets
consistent_junctions <- slope_data_consistent %>%
  filter(Consistent == TRUE) %>%
  select(X)

opposite_junctions <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  select(X)

# Create a contingency table
junctions <- unique(c(consistent_junctions$X, opposite_junctions$X))
consistent_counts <- sapply(junctions, function(j) sum(consistent_junctions$X == j))
opposite_counts <- sapply(junctions, function(j) sum(opposite_junctions$X == j))

contingency_table <- matrix(c(
  sum(consistent_counts > 0), sum(consistent_counts == 0),
  sum(opposite_counts > 0), sum(opposite_counts == 0)
), nrow = 2, byrow = TRUE)

# Perform Chi-Square Test of Independence
chi_square_test_result <- chisq.test(contingency_table)

# Print the results
cat("Chi-Square Test p-value:", chi_square_test_result$p.value, "\n")

# -

# Perform similar operations for the opposite junctions
consistent_p_values <- slope_data_consistent%>%
  filter(Consistent == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now
head(consistent_p_values)

# +

# Perform similar operations for the opposite junctions
opposite_p_values <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now
head(opposite_p_values)

# +
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Calculate summary statistics for consistent test p-values
consistent_p_values <- stat_results_consistent_df %>%
  mutate(Significant = ifelse(mean_p_value < 0.05, "Yes", "No"))

summary_stats_consistent <- consistent_p_values %>%
  summarize(
    mean_p_value_mean = mean(mean_p_value),
    median_p_value_mean = median(mean_p_value),
    sd_p_value_mean = sd(mean_p_value),
    mean_p_value_median = mean(median_p_value),
    median_p_value_median = median(median_p_value),
    sd_p_value_median = sd(median_p_value)
  )

print(summary_stats_consistent)

# +
# Visualize the distribution of mean p-values
ggplot(consistent_p_values, aes(x = mean_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Mean P-values for Consistent Test", x = "Mean P-value", y = "Frequency")

# Visualize the distribution of median p-values
ggplot(consistent_p_values, aes(x = median_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Median P-values for Consistent Test", x = "Median P-value", y = "Frequency")

# QQ plot for mean p-values
ggplot(consistent_p_values, aes(sample = mean_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Mean P-values for Consistent Test")

# QQ plot for median p-values
ggplot(consistent_p_values, aes(sample = median_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Median P-values for Consistent Test")

# Identify significant genes based on mean p-values
significant_genes_mean_consistent <- consistent_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)


# -

length(slope_data_consistent$Consistent)
consistent_table <- table(slope_data_consistent$Consistent)
dim(consistent_table)

# +

# Perform similar operations for the opposite junctions
opposite_p_values <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)

# Create a contingency table for Chi-Square Test
contingency_table <- matrix(c(
  nrow(significant_genes_mean_consistent), nrow(consistent_junctions) - nrow(significant_genes_mean_consistent),
  nrow(significant_genes_mean_opposite), nrow(opposite_junctions) - nrow(significant_genes_mean_opposite)
), nrow = 2, byrow = TRUE)

# Perform Chi-Square Test of Independence
chi_square_test_result <- chisq.test(contingency_table)

# Print the results
cat("Chi-Square Test p-value:", chi_square_test_result$p.value, "\n")


# +
# Visualize the distribution of mean p-values (assuming p-values are available)
consistent_p_values <- slope_data_consistent %>%
  filter(Consistent == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now

# Identify significant genes based on mean p-values
significant_genes_mean_consistent <- consistent_p_values %>%
  filter(mean_p_value < 0.05) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value)

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (consistent test):", nrow(significant_genes_mean_consistent), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_consistent, "significant_genes_mean_consistent.csv", row.names = FALSE)



# +
# Perform similar operations for the opposite junctions
opposite_p_values <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)

# Create a contingency table for Chi-Square Test
contingency_table <- matrix(c(
  nrow(significant_genes_mean_consistent), nrow(consistent_junctions) - nrow(significant_genes_mean_consistent),
  nrow(significant_genes_mean_opposite), nrow(opposite_junctions) - nrow(significant_genes_mean_opposite)
), nrow = 2, byrow = TRUE)

# Perform Chi-Square Test of Independence
chi_square_test_result <- chisq.test(contingency_table)

# Print the results
cat("Chi-Square Test p-value:", chi_square_test_result$p.value, "\n")


# +
# Perform similar operations for the opposite junctions
opposite_p_values <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)

# Create a contingency table for Chi-Square Test
contingency_table <- matrix(c(
  nrow(significant_genes_mean_consistent), nrow(consistent_junctions) - nrow(significant_genes_mean_consistent),
  nrow(significant_genes_mean_opposite), nrow(opposite_junctions) - nrow(significant_genes_mean_opposite)
), nrow = 2, byrow = TRUE)

# Perform Chi-Square Test of Independence
chi_square_test_result <- chisq.test(contingency_table)

# Print the results
cat("Chi-Square Test p-value:", chi_square_test_result$p.value, "\n")


# +
# Identify significant junctions in consistent and opposite sets
consistent_junctions <- slope_data_consistent %>%
  filter(Consistent == TRUE) %>%
  select(X)

opposite_junctions <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  select(X)

# Merge consistent junctions with gene names from featureData
consistent_p_values <- slope_data_consistent %>%
  filter(Consistent == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now


# +

# Visualize the distribution of mean p-values
ggplot(consistent_p_values, aes(x = mean_p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Mean P-values for Consistent Test", x = "Mean P-value", y = "Frequency")

# Visualize the distribution of median p-values
ggplot(consistent_p_values, aes(x = mean_p_value)) + # Assuming mean_p_value for visualization
  geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Histogram of Median P-values for Consistent Test", x = "Median P-value", y = "Frequency")

# QQ plot for mean p-values
ggplot(consistent_p_values, aes(sample = mean_p_value)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Mean P-values for Consistent Test")

# QQ plot for median p-values
ggplot(consistent_p_values, aes(sample = mean_p_value)) + # Assuming mean_p_value for QQ plot
  geom_qq() +
  geom_qq_line() +
  labs(title = "QQ Plot of Median P-values for Consistent Test")

# Identify significant genes based on mean p-values
significant_genes_mean_consistent <- consistent_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Identify significant genes based on median p-values
significant_genes_median_consistent <- consistent_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, median_p_value = mean_p_value) # Assuming same p-value column

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (consistent test):", nrow(significant_genes_mean_consistent), "\n")
cat("Number of significant junctions based on median p-value (consistent test):", nrow(significant_genes_median_consistent), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_consistent, "significant_genes_mean_consistent.csv", row.names = FALSE)
write.csv(significant_genes_median_consistent, "significant_genes_median_consistent.csv", row.names = FALSE)

# Perform similar operations for the opposite junctions
opposite_p_values <- slope_data_opposite %>%
  filter(Opposite == TRUE) %>%
  mutate(gene_name = featureData[X, 2]) %>%
  select(X, gene_name, mean_p_value = Paired_Slope) # Assuming mean_p_value is Paired_Slope for now

# Identify significant genes based on mean p-values
significant_genes_mean_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, mean_p_value)

# Identify significant genes based on median p-values
significant_genes_median_opposite <- opposite_p_values %>%
  filter(mean_p_value < 0.05) %>%
  select(X, gene_name, median_p_value = mean_p_value) # Assuming same p-value column

# Print the number of significant genes
cat("Number of significant junctions based on mean p-value (opposite test):", nrow(significant_genes_mean_opposite), "\n")
cat("Number of significant junctions based on median p-value (opposite test):", nrow(significant_genes_median_opposite), "\n")

# Save significant genes to CSV files
write.csv(significant_genes_mean_opposite, "significant_genes_mean_opposite.csv", row.names = FALSE)
write.csv(significant_genes_median_opposite, "significant_genes_median_opposite.csv", row.names = FALSE)

# Create a contingency table for Chi-Square Test
contingency_table <- matrix(c(
  nrow(significant_genes_mean_consistent), nrow(consistent_junctions) - nrow(significant_genes_mean_consistent),
  nrow(significant_genes_mean_opposite), nrow(opposite_junctions) - nrow(significant_genes_mean_opposite)
), nrow = 2, byrow = TRUE)

# Perform Chi-Square Test of Independence
chi_square_test_result <- chisq.test(contingency_table)

# Print the results
cat("Chi-Square Test p-value:", chi_square_test_result$p.value, "\n")

# -


