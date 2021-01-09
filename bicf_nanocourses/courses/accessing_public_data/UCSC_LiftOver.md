# UCSC Lift Over

[LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)

Must have a bed file
[Bed Format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

Must have at least 3 tab-delimited columns. However, bed 4, bed 5, bed 6, and bed 12 are accepted

Take the ALB.gtf file and convert to bed format. You could use execel to make conversion. But today we will not be using Excel

Please note that the gtf file is 1 base and bed format is 0 base; that means that the gtf sequence starts at 1, and the bed file starts at 0. Both files have the same end point.
```console
[train00@Nucleus000 ~] $ cd ~
[train00@Nucleus000 ~] $ awk '{OFS = "\t"} {print $1, $4-1, $5}' ALB.gtf >ALB.bed
```

On the UCSC Lift Over website, convert human GRCh38 to Human hg19 by choosing "ALB.bed" and submitting the files.

Under Results, select "View Conversions", this downloads a file of hg19 coordinates for ALB transcript

### Workshop
Download the gtf file of ALB in hg19 from UCSC
Compare the coordinates