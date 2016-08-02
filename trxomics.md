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
DESeq2 enables us to use a single function to do this. This function will print out a message for the various steps it performs.
```R
dds <- DESeq(dds)
res <- results(dds)
res
mcols(res, use.names=TRUE)
summary(res)
```
You can also apply some other threshold (smaller p value or bigger logFoldChange value to filter the resutls)
```R
outdf<-cbind(gene_name = rownames(res), data.frame(res))
write.table(outdf,"deseq2.res.xls",quote=F,sep="\t",row.names=F)
head(outdf)

deG <- outdf[(abs(outdf$log2FoldChange)>=1 & outdf$padj<=0.01),]
write.table(deG,"deseq2.deG.xls",quote=F,sep="\t",row.names=F)
dim(deG)

```

Draw heatmap
```R
rld_df <-data.frame(assay(rld))
deG_rld <-rld_df[rownames(rld_df %in% deG$gene_name，]
pheatmap(deG_rld,scale="row",show_rownames = F)
```



###Run StringTie
We will only try to do the transcripts quantification.

In order to use Ballgown program afterwards, "-B" should be used so stringTie will generate Ballgown readable results.

```shell
mkdir stringtie_out
cd stringtie_out
mkdir SRR1551069
mkdir SRR1551068
mkdir SRR1551055
mkdir SRR1551054
mkdir SRR1551048
mkdir SRR1551047
mkdir SRR1550987
mkdir SRR1550986
cd ../
pwd
#You should at your home directory

/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551069.dedup.bam -o stringtie/SRR1551069/SRR1551069.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551068.dedup.bam -o stringtie/SRR1551068/SRR1551068.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551055.dedup.bam -o stringtie/SRR1551055/SRR1551055.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551054.dedup.bam -o stringtie/SRR1551054/SRR1551054.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551048.dedup.bam -o stringtie/SRR1551048/SRR1551048.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1551047.dedup.bam -o stringtie/SRR1551047/SRR1551047.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1550987.dedup.bam -o stringtie/SRR1550987/SRR1550987.gtf
/data/bootcamp/software/stringtie-1.2.4.Linux_x86_64/stringtie -B -G /data/bootcamp/refdb/gencode.gtf -p 6 /data/bootcamp/day3/SRR1550986.dedup.bam -o stringtie/SRR1550986/SRR1550986.gtf


```
Files output in Ballgown readable format are:

+ e_data.ctab: exon-level expression measurements. One row per exon. Columns are e_id (numeric exon id), chr, strand, start, end (genomic location of the exon), and the following expression measurements for each sample:
  + rcount: reads overlapping the exon
  + ucount: uniquely mapped reads overlapping the exon
  + mrcount: multi-map-corrected number of reads overlapping the exon
  + cov average per-base read coverage
  + cov_sd: standard deviation of per-base read coverage
  + mcov: multi-map-corrected average per-base read coverage
  + mcov_sd: standard deviation of multi-map-corrected per-base coverage
+ i_data.ctab: intron- (i.e., junction-) level expression measurements. One row per intron. Columns are i_id (numeric intron id), chr, strand, start, end (genomic location of the intron), and the following expression measurements for each sample:
  + rcount: number of reads supporting the intron
  + ucount: number of uniquely mapped reads supporting the intron
  + mrcount: multi-map-corrected number of reads supporting the intron
+ t_data.ctab: transcript-level expression measurements. One row per transcript. Columns are:
  + t_id: numeric transcript id
  + chr, strand, start, end: genomic location of the transcript
  + t_name: Cufflinks-generated transcript id
  + num_exons: number of exons comprising the transcript
  + length: transcript length, including both exons and introns
  + gene_id: gene the transcript belongs to
  + gene_name: HUGO gene name for the transcript, if known
  + cov: per-base coverage for the transcript (available for each sample)
  + FPKM: Cufflinks-estimated FPKM for the transcript (available for each sample)
+ e2t.ctab: table with two columns, e_id and t_id, denoting which exons belong to which transcripts. These ids match the ids in the e_data and t_data tables.
+ i2t.ctab: table with two columns, i_id and t_id, denoting which introns belong to which transcripts. These ids match the ids in the i_data and t_data tables.


While the program is running, make a folder called isoform, and download the stringtie_out in day3 and "design.pe.txt" to that folder. 
Now let's take a look at the folder structures.
For Ballgown to read the folder, it has to be all the ctab of each sample in an individual sample folder, and all sample folder have similar naming scheme, in our case "SRR155", for Ballgown to identify all the samples.


####Run Ballgown

Install ballgown:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
library("ballgown")
```

Now use the session buttion to change the working directory to "isoform".
```R
stdir <- 'stringtie_out'
design <- "design.pe.txt"

samtbl <- read.table(file=design,header=TRUE,sep='\t')
bg <- ballgown(dataDir=stdir, samplePattern='SRR155', meas='all')
exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
transcript_name <- transcriptNames(bg)
rownames(transcript_gene_table) <- transcript_name
genes <- unique(cbind(geneNames(bg),geneIDs(bg)))
colnames(genes) <- c('SYMBOL','ENSEMBL')

samples <- sampleNames(bg)
mergetbl <- merge(as.data.frame(samples),samtbl,by.x="samples",by.y="File",all.x=TRUE,sort=FALSE)
pData(bg) = data.frame(id=sampleNames(bg), group=as.character(mergetbl$Group),subj=as.character(mergetbl$SubjID))
phenotype_table = pData(bg)
                                        
stat_results = stattest(bg, feature='transcript', meas='cov', covariate='group',adjustvars='subj',getFC=TRUE)

r1 <- merge(transcript_gene_table,na.omit(stat_results),by.y='id',by.x='t_id',all.y=TRUE,all.x=FALSE)
r2 <- merge(genes,r1,by.y='g_id',by.x='ENSEMBL',all.y=TRUE,all.x=FALSE)
write.table(r2,file='de_altsplice.bg.txt',quote=FALSE,sep='\t',row.names=FALSE)


plotMeans('ENSG00000101307.12', bg, groupvar='group', meas='cov', colorby='transcript')



plotLatentTranscripts(gene='ENSG00000101307.12', gown=bg, k=2, method='kmeans', returncluster=FALSE)

```
