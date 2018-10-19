# Workshop from Alignment to Differential Expression

Today we are going to:

* Explore GTF File
* Run featureCounts to get gene counts
* Use edgeR to normalize data and visualization
* Use edgeR to do differential expressed genes
* Differential expressed genes filtering and draw heatmap

## Log into BioHPC
First, we will log into the log into a compute node, install a necessary package and then log into a compute node

* Set up a [WebGui](https://portal.biohpc.swmed.edu/terminal/webgui) session on BioHPC
* Log into via the VNC
* Open a terminal window -- you should be in the directory /archive/nanocourse/June2018/trainXX
* Copy data from /archive/nanocourse/June2018/shared/session2 to your train directory


~~~~
ln -s /archive/nanocourse/June2018/shared/session2 /archive/nanocourse/June2018/trainXX/session2
ls /archive/nanocourse/June2018/trainXX/session2
~~~~

## Exploring GTF File

* How many lines in the GTF file?

* How would you view the first 5 lines of the GTF file?


* What are the first 5 lines of the file?

> ##description: evidence-based annotation of the mouse genome (GRCm38), version M10 (Ensembl 85)
> ##provider: GENCODE
> ##contact: gencode-help@sanger.ac.uk
> ##format: gtf
> ##date: 2016-07-19

* Why might this information be important?

> Indicates the version of the annotation use and when it was generated.

* How many transcripts are there for Tnnt2?
(hint: use grep, cut, uniq 9th field (attributes) transcript ids
have ENSMUST in transcript_id field semicolon ’;’ can be used as a delimiter)


* How many transcripts are protein coding (hint: transcript_type)?

## Run FeatureCounts to get counts on a file

```
module load subread/1.6.1
featureCounts -p -g gene_name -a session2/gencode.gtf -o heart_e11_5_rep1.cts session2/alignments/heart_e11_5_rep1.dedup.bam
```

featureCounts will print statistics on the screen. There is one line in the second box:
"Successfully assigned fragments : 376082 (72.3%)"

Depending on how your samples are prepared, this value varies:
* If you are using poly-A extraction, we expect this value to be at least 70%
* If you are extracting total RNA, this value will be around 50%
* When you have a percent less than 40%, you should look into this matter to figure out why

Now Let's take a look at featureCounts results
```shell
head heart_e11_5_rep1.cts.summary
head heart_e11_5_rep1.cts
```

* What are the counts for the following genes?
  * Tnnt2
  * Myh6
  * Myh7
  * Kcng2  
  * Gja5  

Run feature counts on all of the files.


Then we need to combine all the counts into a matrix.
We will be using an in-house Perl script "concat_cts.pl" to do this.
It takes in the names of all genes (genenames.txt)

There are 2 parameter you need to provide:
1. The output folder -o (here we use the current directory ./ )
2. the patter of all of the cts files (.cts)

```
cp session2/genenames.txt .
perl session2/scripts/concat_cts.pl -o ./ session2/counts/*.cts
head countTable.txt
```

Before we are ready for differential expression we need to create
a design file. Thank full we have already done that for you.

Also lets copy 2 files to your home directory

```
cp session2/design_se.txt /home2/trainXX
cp countTable.txt /home2/trainXX
```



## Use RStudio on demand

* Set up a [RStudio](http://rstudio.biohpc.swmed.edu/) session on BioHPC


### Install edgeR for gene differential expression analysis
Install edgeR

```R
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```
For pheatmap package, click "package"->"install" Type 'a' if program asks you if you want to update packages.
Install from:CRAN Packages: pheatmap
Check "Install dependencies"
Click "Install"

Click "Yes" if prompt window asks you if you want to use a personal library.

After everything finished, load the libraries using
```R
library("edgeR")
library("pheatmap")
```

### Load Data


```R

#Read data matrix and sample file

cfile<-read.table("countTable.txt",header=T,row.names=1)
coldata<-read.table("design_se.txt",header=T,sep="\t")

head(counts)
head(coldata)

#Reorder the counts columns to match the order of sample file
cfile = cfile[c("heart_e11_5_rep1", "heart_e11_5_rep2", "heart_p0_rep1", "heart_p0_rep2")]

#It is good to set your control group label as the baseline. Especially you are going to use intercept
group = relevel(factor(coldata$SampleGroup),ref="heart_e11_5")
cds = DGEList(cfile,group=group)

```
### Pre-filtering the low-expressed genes

Filter for keeping a gene if cpm (counts per million) exceeds 1 in at least 2 samples.
```
cds = cds[ rowSums(cpm(cds)>=1) >= 2, ,keep.lib.sizes=FALSE]
```

#### Exploratory analysis and some vizulization

Use cpm() function to get log2 transformed normalized counts

```
rld <- cpm(cds, log=TRUE)
```

Calculate the distance between sample pairs and do hierarchical clustering
```
sampleDists = dist(t(rld))
sampleDists
plot(hclust(sampleDists))
```


Use heatmap to show sample correlation**

```
pheatmap(cor(rld))
```

Use MDS plot to check the relationship of replicates.
```
points = c(15,16)
colors = rep(c("red","blue"),4)
plotMDS(cds, col=colors[group], pch=points[group])
legend("bottomleft", legend=levels(group), pch=points, col=colors, ncol=2)
```

Make a design matrix for samplegroups
```
samplegroup <- factor(coldata$SampleGroup)
design<-model.matrix(~samplegroup)
design
```

Normalize data and estimate dispersion. What is the norm.factors per sample?

```
cds = calcNormFactors(cds)
cds
```

Use glmFit() and glmLRT() to test for differential expression. What are the top
10 differential genes sorted by logFC?

```
cds$samples
cds <- estimateDisp(cds,design)

fit <- glmFit(cds, design)
lrt <- glmLRT(fit, coef=2)
res <- topTags(lrt, n=dim(cfile)[1],sort.by="logFC")
res[1:10,]
```

Make a dataframe with column for genes

```
res_df = cbind(gene_name = rownames(res), data.frame(res))
write.table(res_df,"edgeR.res.tsv",quote=F,sep="\t",row.names=F)
head(res_df)
```

Filter for genes that are logFC >= 1 & FDR <= 0.01
```
res_filt = res_df[(abs(res_df$logFC)>=1 & res_df$FDR<=0.01),]
```

Draw a heatmap of differential expressed genes meeting filter criteria.
```
res_filt_rld = rld[rownames(rld) %in% res_filt$gene_name,]
pheatmap(res_filt_rld,scale="row",show_rownames = F)
```
