# SRA Downloads

[SRA webpage](https://www.ncbi.nlm.nih.gov/sra)

## Option1: Download SRA data directly from website



#### Get a list of SRA ID's to download
Find the data in SRA you would like

Search: 
```
lymph node tissue in BALB/c mice
```

Click on replicate 3. There are 8 SRR files, single-end, RNA-seq.

In the upper right-hand corner,<br/>
&nbsp;&nbsp;select "Send to",<br/>
&nbsp;&nbsp;select "File",<br/>
&nbsp;&nbsp;under format select "Accession List",<br/>
&nbsp;&nbsp;then finally select "Create File"

Now you have a text file list of SRA ID's to download.



#### Use SRA Toolkit to download SRA files

###### Download the SRA program (not necessary for this class)
[SRA Toolkit Software](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software), 
Supports Linux, Unix, and Mac
(Don't do this here in class, we will use BioHPC's install)

###### Run SRA program on BioHPC
Open up a terminal in your webvizualization
cd into your working folder /work/archive/nanocourse/access_pub_data/train00
```console
[train00@Nucleus000 ~] $ cd ~
```

Make and enter a folder for SRA
```console
[train00@Nucleus000 ~] $ mkdir SRA_training
[train00@Nucleus000 ~] $ cd SRA_training
```


Load the module for SRA (only necessary when on a cluster)
```console
[train00@Nucleus000 ~] $ module load sra_toolkit/2.8.2-1
```

Find the path to your file
```console
[train00@Nucleus000 ~] $ srapath SRR9729388
https://sra-download.ncbi.nlm.nih.gov/traces/sra5/SRR/009501/SRR9729388
```

Use wget to download file. This is not the same as the manual. Due to experience not all SRA files download properly with SRA prefetch
```console
[train00@Nucleus000 ~] $ wget https://sra-download.ncbi.nlm.nih.gov/traces/sra5/SRR/009501/SRR9729388 
```

Confirm that the download was correct and complete
```console
[train00@Nucleus000 ~] $ vdb-validate SRR9729388
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Database 'SRR9729388' metadata: md5 ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Table 'SEQUENCE' metadata: md5 ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'ALTREAD': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'QUALITY': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'READ': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'READ_LEN': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'READ_START': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Column 'SPOT_GROUP': checksums ok
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Database 'SRR9729388' contains only unaligned reads
2020-01-21T20:21:59 vdb-validate.2.8.2 info: Database 'SRR9729388' is consistent
```


Use fastq-dump to transform SRA file into fastq files<br/>
Options:<br/>
&nbsp;&nbsp;--gzip: Compress output using gzip<br/>
&nbsp;&nbsp;--split-3: Legacy 3-file splitting for mate-pairs: First biological reads satisfying dumping conditions are placed in files *_1.fastq and *_2.fastq If only one biological read is present it is placed in *.fastq Biological reads and above are ignored.<br/>
```console
[train00@Nucleus000 ~] $ fastq-dump --gzip --split-3 `readlink -e SRR9729388`
```
If you do not use the "readlink -e " in the command, a SRA cache will be created somewhere on your computer taking up valueable resources.







## Option2: Use Astrocyte to download files for you (requires BioHPC account)

Make a design file
"This pipeline requires a very simple design, tab delimited file: Line 1 must be the following header: "sample_id sra_number" Each subsequent line should then contain, in order the sample ID (whatever you wish to call the sample), and the SRA number from NCBI."
One has been created for you and in /work/archive/nanocourse/access_pub_data/shared


[BioHPC website](https://portal.biohpc.swmed.edu)

Under "Cloud Services", click "Astrocyte Workflow Platform"

In the upper right hand corner, login to your training account.

Select "Start a Project"

At the top, create a project name. For instance "SRA_training". Then click "Create New Project".

Click "Add Data To This Project"

Under "Upload files from the web", click "Select file to upload..."

Navagate to the course shared folder and select "SRA_design.tsv"

Click "Finished uploading files"

Under a header "Worflows run in this project", click "Run a workflow in this project".

Select "Run Workflow" next to the "BICF SRA Download Pipeline"

Set up the run by first making sure that the Project and design files are correct. Then add a "Name for this run"; for instance "SRAtest".

At this time you can close the web browser without losing any information. 






## Workshop
Find, download, validate download, and convert SRA to fastq<br/>
1) The species or cell line of your current project<br/>
2) RNA-seq or ChiP-seq<br/>
3) Optional: Organ, Tissue or Anatomic Site<br/>

