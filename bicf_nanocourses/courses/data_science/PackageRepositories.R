# Load the count data

# Create group labels
# Rownames = sample names
# Groups column with labels
# NoBCH = S085BE, S089BE, S114BE, S118BE
# BCH = S110BE, S117BE, S132INBE


# Data cleaning
# Filter data to only "protein_coding" genes
# Filter to only the firs 1000 genes (to speed up this workshop)
# Set the rownames as the gene symbols
# Remove the unnecessary columns (1, 2, and 3)
# Convert to a matrix


# Plot a heatmap with base R

# Replot with no rownames

# Replot scaled by gene (make a separate matrix because we will use it later)


# Install CRAN package pheatmap

# Use pheatmap to plot scaled data

# Replot without row names


# Install ComplexHeatmap from the GitHub repository https://github.com/jokergoo/ComplexHeatmap

# Use ComplexHeatmap to plot scaled data

# Replot without row names


# Install DESeq2 from Bioconductor

# Load the data into the DESeqDataSet class

# Run DESeq2

# Extract DESeq results into DESeqResults class

# Extract a list of genes which are significant (alpha = 0.05)


# Draw a heatmap of only the significan genes
