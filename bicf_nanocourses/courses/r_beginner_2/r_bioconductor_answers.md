# R Bioconductor RNA-seq and ChIP-seq analysis

## General Software questions

In this exercise, we are going to download two packages and explore
some of their functionality.

1. How do you download any package from Bioconductor ?
```
source("https://bioconductor.org/biocLite.R")
biocLite("PackageName")
```
2. How do you access the vignettes for a given package ?
```
browseVignettes("PackageName")
```

Download [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
and [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html).

3. What is the current version of ChIPseeker?
```
1.14.1 (found on ChIPseeker Page)
```
4. How should you cite edgeR?
```
Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), pp. 139-140.

McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” Nucleic Acids Research, 40(10), pp. 4288-4297.
```
5. How do you find help for a function (say annotatePeak)?
```
?annotatePeak
```

## ChIP-seq analysis using ChIPseeker

In this exercise, we are going to explore annotating, comparing
and visualizing ChIP-seq data and working with GRanges objects.
Open the vignette for ChIPSeeker to load test data.

6. How many peaks are in the ARmo_100nM experiment?
```
peak_ARmo_100nM <- readPeakFile(files$ARmo_100nM)
length(peak_ARmo_100nM)
```
7. What is the average length of peaks in ARmo_1nM?
```
peak_ARmo_1nM <- readPeakFile(files$ARmo_1nM)
mean(width(peak_ARmo_1nM))
```
8. Extract peaks only from chromosome 3 in ARmo_1nM?
```
peak_ARmo_1nM_chr3 <- peak_ARmo_1nM[seqnames(peak_ARmo_1nM) == 'chr3']
```
9. Extract the first five peaks from the ARmo_1nM?
```
peak_ARmo_1nM[1:5]
```
10. Extract the peak name and score column of CBX6_BF?
```
peak_CBX6_BF <- readPeakFile(files$CBX6_BF)
mcols(peak_CBX6_BF)
```
11. What is the most represented annotation of peaks in CBX7_BF?
```
peakAnno_CBX7_BF <- annotatePeak(files$CBX7_BF, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_CBX7_BF
```
12. Make a vendiagram of overlapping peaks for ARmo experiments?
```
ARmo_files <- files[1:3]
ARmo_labels = c('ARmo_0M', 'ARmo_1nM', 'ARmo_100nM')
vennplot.peakfile(ARmo_files, labels = ARmo_labels)
```
13. Get the Average profile around promoters for CBX7_BF?
```
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(files$CBX7_BF, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000))
```

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
```
cds = cds[ rowSums(cpm(cds)>=1) >= 4, ,keep.lib.sizes=FALSE]
```

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
```
points = c(15,16)
colors = rep(c("red","blue"),4)
plotMDS(cds, col=colors[group], pch=points[group])
legend("bottomleft", legend=levels(group), pch=points, col=colors, ncol=2)
```
16. Make a design matrix for samplegroups
```
samplegroup <- factor(coldata$SampleGroup)
design<-model.matrix(~samplegroup)
design
```
17. Normalize data and estimate dispersion. What is the norm.factors per sample?
```
cds = calcNormFactors(cds)
cds
```
18. Use glmFit() and glmLRT() to test for differential expression. What are the top
10 differential genes sorted by logFC.
```
cds$samples
cds <- estimateDisp(cds,design)

fit <- glmFit(cds, design)
lrt <- glmLRT(fit, coef=2)
res <- topTags(lrt, n=dim(cfile)[1],sort.by="logFC")
res[1:10,]
```

**Make a dataframe with column for genes**
```
res_df = cbind(gene_name = rownames(res), data.frame(res))
```

19. Filter for genes that are logFC>=1 & FDR<=0.01
```
res_filt = res_df[(abs(res_df$logFC)>=1 & res_df$FDR<=0.01),]
```
20. Draw a heatmap of differential expressed genes meeting filter criteria.
```
res_filt_rld = rld[rownames(rld) %in% res_filt$gene_name,]
pheatmap(res_filt_rld,scale="row",show_rownames = F)
```
