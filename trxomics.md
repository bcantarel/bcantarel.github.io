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
```bash
mkdir day3
cp /data/bootcamp/day3/* day3/
cd day3
ls
```
And you should see:
```

design.pe.txt         SRR1551047.dedup.bam  SRR1551055.dedup.bam
SRR1550986.bam        SRR1551048.bam        SRR1551068.bam
SRR1550986.cts        SRR1551048.cts        SRR1551068.cts
SRR1550986.dedup.bam  SRR1551048.dedup.bam  SRR1551068.dedup.bam
SRR1550987.bam        SRR1551054.bam        SRR1551069.bam
SRR1550987.cts        SRR1551054.cts        SRR1551069.cts
SRR1550987.dedup.bam  SRR1551054.dedup.bam  SRR1551069.dedup.bam
SRR1551047.bam        SRR1551055.bam
SRR1551047.cts        SRR1551055.cts

```

Then: remove the PCR duplicates in your samples (We will only try one sample)
```sh
/data/bootcamp/seqprg/samtools-1.3.1/samtools rmdup day3/SRR1551047.bam  test.bam
```

Type: 
```sh
/data/bootcamp/software/subread-1.5.0-p3-source/bin/featureCounts
```
You should get help message of featureCounts on your screen.

