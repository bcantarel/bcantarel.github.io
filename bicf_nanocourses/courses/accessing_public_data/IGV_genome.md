# Make a custom genomve for IGV web app

[IGV web app](https://igv.org/app/)


### Load the fasta

In IGV select, "Genome", "Local File...", "ALB.fa"

```console
Genome File Error

The following files were not loaded...
ALB.fa ERROR: index file must also be selected
```

To load a custom genome, you must make an index file
Make an index file
```console
[train00@Nucleus000 ~] $ cd ~
[train00@Nucleus000 ~] $ module load samtools
[train00@Nucleus000 ~] $ samtools faidx ALB.fa
```

Now load the fasta file and index, together, into IGV


### Load the gtf file into IGV

In IGV select, "Tracks", "Local File...", "ALB.gtf"

Notice that there isn't an error, but the screen is blank. We would expect if the file properly loaded that there would be colored blocks of annotations on the screen. This is because we need to make the fasta and gtf match.

Here are the mismatches:
1) The fasta calls the chromosome name "hg38_ncbiRefSeq_NM_000477.7" and the gtf calls it "chr4"
2) The sequence isn't annotated for start to end so it starts at base 1 and ends at base 19,246. Whereas, the gtf file annotates the chromosomal locations ranging from 73,404,239 to 73421484. The bases need to match.

Fix the gtf file:
DO NOT OPEN THE GTF FILE IN A WORD DOCUMENT!

On your webviz session
&nbsp;&nbsp;select "Applications", "Accessories", "Text Editor"
&nbsp;&nbsp;select "Open", "ALB.gtf"
&nbsp;&nbsp;In the file, replace "chr4" to "hg38_ncbiRefSeq_NM_000477.6"
&nbsp;&nbsp;In the file, subtract "73404238" (one less than the smallest number) on columns 4 and columns 5

This would be easier in excel but we don't have it. Here's the command to fix it
```console
[train00@Nucleus000 ~] $ awk '{OFS = "\t"} {print "hg38_ncbiRefSeq_NM_000477.7", $2, $3, $4-73403238, $5-73403238, $6, $7, $8, $9}' ALB.gtf >ALB2.gtf
```


### Remove the bad gtf files and save image for publication

On the right of each track is settings, click the wheel, click "Remove track"
To save, click "Save SVG", "OK"
