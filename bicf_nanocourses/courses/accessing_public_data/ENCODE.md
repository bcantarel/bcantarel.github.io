# ENCODE Downloads

[ENCODE Project](https://www.encodeproject.org/)


### ENCODE Homepage
Start your search

### Data: Experiment Matrix
Refine your search

### Encyclopedia: Visualize (SCREEN)
SCREEN: Search Candidate cis-Regulatory Elements by ENCODE

### Materials & Methods: Antibodies
Confirm lot number and efficacy of antibodies

### Materials & Methods: Assays and standards
Description of standard pipeline.
Minimal quality control metrics for experiments.




#### Download data, single item
On the website, find the data for ENCFF002CTW

Create and enter a new folder
```console
[train00@Nucleus000 ~] $ cd ~
[train00@Nucleus000 ~] $ mkdir encode_test
[train00@Nucleus000 ~] $ cd encode_test
```

Single project
```console
[train00@Nucleus000 ~] $ wget https://www.encodeproject.org/files/ENCFF002CTW/@@download/ENCFF002CTW.bed.gz
```

Validate download is complete with MD5sum and content MD5sum
```console
[train00@Nucleus000 ~] $ md5sum ENCFF002CTW.bed.gz 
92ef073fbaf203c55bf9a8420c1d47e8  ENCFF002CTW.bed.gz

[train00@Nucleus000 ~] $ zcat ENCFF002CTW.bed.gz | md5sum
d1f0e92a2b94a7160c1a4ca90a54c560  -
```

#### Download data, bulk item/shopping cart
On the website, find the data for Jurkat, ChiP-seq, replicated peaks


From shopping cart<br/>
&nbsp;&nbsp;download files.txt (open in text editor and then save in encode_test folder)
&nbsp;&nbsp;do a batch download of all files in list
```console
[train00@Nucleus000 ~] $ xargs -L 1 curl -O -L < files.txt
```

The files and their metadata are saved as a output
Get the md5sum from the metadata (can open metadata in excel)
```console
[train00@Nucleus000 ~] $ awk -F "\t" '{OFS = "\t"} {print $1, $41}' \?type\=Experiment\&files.output_type\=replicated+peaks
```

Check the MD5sum of each file
```console
[train00@Nucleus000 ~] $ md5sum *
```

The content md5sum is NOT found in the meta data
Optional: Check the content MD5sum (requires loop)
```console
[train00@Nucleus000 ~] $ for i in `ls`; do zcat ${i} | md5sum; done
```
Notice the error in the output


## Workshop 1
Find and download the fastq file of
1) The species or cell line of your current project<br/>
2) RNA-seq or ChiP-seq<br/>
3) Optional: Organ, Tissue or Anatomic Site<br/>

## Workshop 2
Find a cis-regularoty element in
1) The species or cell line of your current project<br/>
2) Look at the 1st protein-coding gene in the list<br/>

How many other accessions have more expression of that gene than your accession?<br/>
Look at the UCSC browser of up to 5 of those accessions. What are the differences in the pileup of the histone marks?<br/>

