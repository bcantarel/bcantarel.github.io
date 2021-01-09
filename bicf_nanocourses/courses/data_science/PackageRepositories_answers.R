# Load the count data
countTable <- read.delim("countTable.txt")

# Create group labels
# Rownames = sample names
# Groups column with labels
# NoBCH = S085BE, S089BE, S114BE, S118BE
# BCH = S110BE, S117BE, S132INBE
groups <- data.frame(Groups=c("NoBCH","NoBCH","BCH","NoBCH","BCH","NoBCH","BCH"),row.names=colnames(countTable[,4:10]))


# Data cleaning
# Filter data to only "protein_coding" genes
# Filter to only the firs 1000 genes (to speed up this workshop)
# Set the rownames as the gene symbols
# Remove the unnecessary columns (1, 2, and 3)
# Convert to a matrix
data <- countTable
data <- data[data$TYPE=="protein_coding",]
data <- data[1:1000,]
rownames(data) <- data$SYMBOL
data <- data[,-c(1:3)]
data <- as.matrix(data)


# Plot a heatmap with base R
heatmap(data)

# Replot with no rownames
heatmap(data,labRow=FALSE)

# Replot scaled by gene (make a separate matrix because we will use it later)
data_scaled <- t(scale(t(data)))
heatmap(data_scaled,labRow=FALSE)


# Install CRAN package pheatmap
install.packages("pheatmap")
library(pheatmap)

# Use pheatmap to plot scaled data
pheatmap(data_scaled)

# Replot without row names
pheatmap(data_scaled,show_rownames=FALSE)


# Install ComplexHeatmap from the GitHub repository https://github.com/jokergoo/ComplexHeatmap
install.packages("devtools")
#Sys.setenv(http_proxy = "http://proxy.swmed.edu:3128")
#Sys.setenv(https_proxy = "http://proxy.swmed.edu:3128")
devtools::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# Use ComplexHeatmap to plot scaled data
Heatmap(data_scaled)

# Replot without row names
Heatmap(data_scaled,show_row_names=FALSE)


# Install DESeq2 from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load the data into the DESeqDataSet class
dds <- DESeqDataSetFromMatrix(countData=data,colData=groups,design=~Groups)

# Run DESeq2
dds <- DESeq(dds)

# Extract DESeq results into DESeqResults class
res <- results(dds)

# Extract a list of genes which are significant (alpha = 0.05)
res <- res[!is.na(res$padj),]
res <- res[order(res$padj),]
genes_sig <- rownames(res[res$padj <= 0.05,])


# Draw a heatmap of only the significan genes
Heatmap(data_scaled[rownames(data_scaled) %in% genes_sig,])
