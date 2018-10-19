# R Bioconductor RNA-seq and ChIP-seq analysis

## General Software questions

In this exercise, we are going to download two packages and explore
some of their functionality.

1. How do you download any package from Bioconductor ?
2. How do you access the vignettes for a given package ?

Download [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
and [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html).

3. What is the current version of ChIPseeker?
4. How should you cite edgeR?
5. How do you find help for a function (say annotatePeak)?

## ChIP-seq analysis using ChIPseeker

In this exercise, we are going to explore annotating, comparing
and visualizing ChIP-seq data and working with GRanges objects.
Open the vignette for ChIPSeeker to load test data.

6. How many peaks are in the ARmo_100nM experiment?
7. What is the average length of peaks in ARmo_1nM?
8. Extract peaks only from chromosome 3 in ARmo_1nM?
9. Extract the first five peaks from the ARmo_1nM?
10. Extract the peak name and score column of CBX6_BF?
11. What is the most represented annotation of peaks in CBX7_BF?
12. Make a vendiagram of overlapping peaks for ARmo experiments?
13. Get the Average profile around promoters for CBX7_BF?

## RNA-seq analysis using edgeR

In this exercise, we are going to use edgeR to normalize data,
do differential expressed genes, draw heatmaps.
Open the manual and vignette for edgeR for help.

**Install Packages**
For pheatmap package, click "package"->"install" Type 'a' if program asks you if you want to update packages.
Install from:CRAN Packages: pheatmap
Check "Install dependencies"
Click "Install"
Click "Yes" if prompt window asks you if you want to use a personal library.

**Load Data**
```
#Set working directory
#Session-> Set working directory -> choose directory

#Read data matrix and sample file
cfile = read.table("RNAseq_count_table.txt",header=T,row.names=1,sep="\t")
coldata = read.table("RNAseq_design_pe.txt",header=T,sep="\t")
head(cfile)
head(coldata)

#Reorder the counts columns to match the order of sample file
cfile = cfile[coldata$SampleID]

#It is good to set your control group label as the baseline. Especially you are going to use intercept
group = relevel(factor(coldata$SampleGroup),ref="monocytes")
cds = DGEList(cfile,group=group)
```

**Pre-filtering the low-expressed genes.**

14. Filter for keeping a gene if cpm (counts per million) exceeds 1 in at least 4 samples.


#### Explortory analysis and some visulization

Use cpm() function to get log2 transformed normalized counts
```
rld <- cpm(cds, log=TRUE)
```

**Calculate the distance between sample pairs and do hierarchical clustering**
```
sampleDists = dist(t(rld))
sampleDists
plot(hclust(sampleDists))
```


**Use a heatmap to show sample correlation**

```
pheatmap(cor(rld))
```

15. Use MDS plot to check the relationship of replicates.
16. Make a design matrix for samplegroups
17. Normalize data and estimate dispersion. What is the norm.factors per sample?
18. Use glmFit() and glmLRT() to test for differential expression. What are the top
10 differential genes sorted by logFC.

**Make a dataframe with column for genes**
```
res_filt = cbind(gene_name = rownames(res), data.frame(res))
```

19. Filter for genes that are logFC>=1 & FDR<=0.01
20. Draw a heatmap of differential expressed genes meeting filter criteria.
