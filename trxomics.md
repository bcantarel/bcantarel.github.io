##2016 NGS Summer BootCamp: From Alignment to differential expression analysis

####Today we are going to do:
+ Run featureCounts to get gene counts
+ Use DESeq2 to normalize data and visualization
+ Use DESeq2 to do differential expressed genes
+ Differential expressed genes filtering and draw heatmap
+ Run stringTie to get transcript-level abundance
+ Use Ballgown to visualize the transcripts and do differential analysis 


First, we will log into the log into a compute node,install a necessary package and then log into a compute node

Log into your machine or account.

Mac Users: Open the application: Terminal
From Terminal: ssh -Y username@toxea.swmed.edu
Windows users -- Please refer to the PuTTY instructions with your username and the server toxea.swmed.edu


####Run featureCounts
First: copy all the data files we will need today to your home directory
```shell
mkdir day3
cp /data/bootcamp/day3/*cts day3/
cp /data/bootcamp/day3/SRR1551047*bam day3/
cp /data/bootcamp/day3/design.pe.txt day3/
cd day3
ls
```

Optional: remove the PCR duplicates in your samples (This takes a long time to run one sample, please try it by yourself)
```shell
/data/bootcamp/seqprg/samtools-1.3.1/samtools sort -o SRR1551047.sort.bam SRR1551047.bam
/data/bootcamp/seqprg/samtools-1.3.1/samtools rmdup SRR1551047.sort.bam  SRR1551047.dedup.bam
```

Second: try to use featureCounts to get counts 
```shell
/data/bootcamp/software/subread-1.5.0-p3-source/bin/featureCounts
/data/bootcamp/software/subread-1.5.0-p3-source/bin/featureCounts -T 4 -p -g gene_name -a /data/bootcamp/refdb/gencode.gtf -o test.cts day3/SRR1551047.dedup.bam
```

featureCoutns will print statistics on the screen. There is one line in the second box:
"Successfully assigned fragments : 376082 (72.3%)"

Depending on how your samples are prepared, this value varies:
+ If you are using poly-A extraction, we expect this value to be at least 70%
+ If you are extracting total RNA, this value will be around 50%
+ When you have a percent less than 40%, you should look into this matter to figure out why

Now Let's take a look at featureCounts results
```shell
head test.cts.summary
head test.cts
```

Then we need to combine all the counts into a matrix.
We will be using an in-house Python script "combineFeatureCounts.py" to do this.
There are three parameter you need to provide:
1. Input file names. A list contain all the names of the count files;
2. The names of all genes. Here we are going to use one count file to provide this;
3. The name of the output file, which is your count matrix.

```python
ls *.cts > files
python combineFeatureCounts.py files SRR1551047.cts count.table
head count.table
```

Before we start to do the differential expression analysis, please creat a folder "DE_test" on your Desktop or any place you like and use FTP tool to drag "count.table" and "design.txt" into this folder.


###Run DESeq2 for gene differential expression analysis
Open your Rstudio and install DESeq2
```R
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```
For pheatmap package, click "package"->"install"
Type 'a' if program asks you if you want to update packages.
Install from:CRAN
Packages: pheatmap
Check "Install dependencies"
Click "Install"

Click "Yes" if prompt window asks you if you want to use a personal library.

After everything finished, load the libraries using
```R
library("DESeq2")
library("pheatmap")
```

####Load Data
```R
#Read data matrix and sample file
cfile<-read.table("counts.table",header=T,row.names=1)
coldata<-read.table("design.pe.txt",header=T)
head(counts)
head(coldata)

#Reorder the counts columns to match the order of sample file
cfile <- cfile[rownames(coldata)]
dds<-DESeqDataSetFromMatrix(countData=cfile,colData=coldata, design=~Tissue)
```
There are some other things we need to do before the differential analysis.
1. Always check the levels to make sure the control is the first level
```R
dds$Tissue
```
You should see:
"Levels: monocytes neutrophils"
To be safe, it is always good to re-level the factor level:
```R
dds$Tissue <- relevel(dds$Tissue, ref="monocytes")
```

2. Pre-filtering the low-expressed genes. The code below means keep a gene if its counts exceed 10 in at least 4 samples

```R
nrow(dds)
dds <- dds[ rowSums(counts(dds)>=10) >= 4, ]
nrow(dds)
```

####Explortory analysis and some visulization
Use rlog() function to get log2 transformed normalized counts
```R
rld <- rlog(dds, blind=FALSE)
plot(assay(rld)[,1:2],pch=16, cex=0.3)
```

**Calculate the distance between sample pairs and do hierarchical clustering**
```R
sampleDists <- dist( t( assay(rld) ) )
sampleDists
plot(hclust(sampleDists))

```

**PCA plot**
```R
plotPCA(rld, intgroup = c("Tissue"))
```

####Find differential expressed genes
You can use a single function to do this. This function will print out a message for the various steps it performs.
```R
dds <- DESeq(dds)
res <- results(dds)
res
mcols(res, use.names=TRUE)
summary(res)
```
You can also apply some other threshold (smaller p value or bigger logFoldChange value to filter the resutls)
```R
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)

resLFC1 <- results(dds, lfcThreshold=1)
```

For me, I usually save the whole table and do filtering afterwards. 
```R
outdf<-cbind(gene_name = rownames(res), data.frame(res))
write.table(outdf,"deseq2.res.xls",quote=F,sep="\t",row.names=F)
```

###Run StringTie

Install ballgown:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
```

