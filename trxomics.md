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
Type: (Use Tab key to auto fill the address)
```shell
/data/bootcamp/software/subread-1.5.0-p3-source/bin/featureCounts
```
You should get help message of featureCounts on your screen.

