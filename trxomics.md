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


Install ballgown:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
```

