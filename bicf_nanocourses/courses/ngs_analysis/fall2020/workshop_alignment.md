# 2016 NGS Summer BootCamp: Workshop Short Read Alignment

Today we are going to:

* Run FastQC to determine sequence quality
* Trim Sequences
* Align Genomic Sequence to a reference genome
* Align RNA Sequence to a reference genome
* Examine Alignments in IGV

First, we will log into the log into a compute node,install a necessary package and then log into a compute node

* Log into your machine or account.
  1. Mac Users: Open the application: Terminal
  2. From Terminal: `ssh -Y  student01@nucleus.biohpc.swmed.edu
 are going to unzip this file and run fastqc

`
gunzip NA12878_100ng.R*.1M.fq.gz
#FastQC needs an updated version of java
export PATH=/data/bootcamp/seqprg/jre1.8.0_101/bin:$PATH
/data/bootcamp/seqprg/FastQC/fastqc NA12878_100ng.R1.1M.fq
/data/bootcamp/seqprg/FastQC/fastqc  NA12878_100ng.R2.1M.fq
`
Next we will trim the sequences and run fastqc again
`
/data/bootcamp/seqprg/trim_galore_zip/trim_galore --paired -q 25 --illumina --length 35 --no_report_file --path_to_cutadapt  /data/bootcamp/seqprg/bin/cutadapt  NA12878_100ng.R1.1M.fq NA12878_100ng.R2.1M.fq
/data/bootcamp/seqprg/FastQC/fastqc NA12878_100ng.R1.1M_val_1.fq
/data/bootcamp/seqprg/FastQC/fastqc NA12878_100ng.R2.1M_val_2.fq
`
Now we are going to run our first alignment using only the sequences against chrX (for time sake) against GRCh37
*Note: The alignment step (bwa) might take a little while, so once the alignment is running we will skip to the next step, which will entail examining the sequence data
`
ln -s /data/bootcamp/day1/NA12878_100ng.R*.chrX.fastq.gz .
/data/bootcamp/seqprg/bwa-0.7.15/bwa mem -M -R '@RG\tLB:tx\tPL:illumina\tID:NA12878\tPU:barcode\tSM:NA12878' /data/bootcamp/refdb/genome.fa NA12878_100ng.R1.chrX.fastq.gz NA12878_100ng.R2.chrX.fastq.gz > NA12878.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools view -b -u -S -o NA12878.unsort.bam NA12878.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools sort -o NA12878.bam NA12878.unsort.bam
/data/bootcamp/seqprg/samtools-1.3.1/samtools index NA12878.bam
`
While the bwa is aligning the sequences, we will examine the sequence files
If you are using windows, use WINSCP or your favorite data transfer tool to download the fastqc output onto your computer.
If you are using a mac, open a new terminal window and follow these instruction
`
mkdir bootcamp_day1
cd bootcamp_day1
scp username@toxea.swmed.edu:~/day1/\*fastqc\* .
`
Open the html files in your favorite web browser (they should open on double click) and examine all of the quality metrics

First and foremost, the FastQC "Summary" should generally be ignored. Its "grading scale" (green - good, yellow - warning, red - failed) incorporates assumptions for a particular kind of experiment, and is not applicable to most real-world data. Instead, look through the individual reports and evaluate them according to your experiment type.

The FastQC reports I find most useful are:

* The Per base sequence quality report in the raw, which can help you decide if sequence trimming is needed before alignment.
* The Per base sequence quality report in the trimmed, which will show you the affect of trimming
* The Sequence Duplication Levels report, which helps you evaluate library enrichment / complexity.
* The Overrepresented Sequences report, which helps evaluate adapter contamination.

You can have a look through the QC results to try to answer the following questions:

* Did any of the QC modules trigger a warning or alert condition?
* How do the base call qualities provided by the sequencer suggest the data is high quality, or might it benefit from being quality trimmed compare between the raw and trimmed sequence data?
* Are there any consistent sequence biases in the data?
*  Is there any suggestion of the presence of adapter sequence which might need to be removed in the raw data?
*  Does the duplication level of the data look reasonable? 

Once the bwa alignment is complete, we will now do a splice aware alignment, using hisat2 to GRCh38
There are 2 samples, run alignment in both samples: SRR1551068 and SRR1551069
During the alignment step, we will examine the DNA alignment

`
ln -s /data/bootcamp/day1/SRR155106*.fastq.gz .
/data/bootcamp/seqprg/hisat2-2.0.4/hisat2 -p 4 --no-unal --dta -x /data/bootcamp/refdb/genome -1 SRR1551068_1.fastq.gz -2 SRR1551068_2.fastq.gz -S SRR1551068.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools view -b -u -S -o SRR1551068.unsort.bam SRR1551068.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools sort -o SRR1551068.bam SRR1551068.unsort.bam
/data/bootcamp/seqprg/samtools-1.3.1/samtools index SRR1551068.bam

/data/bootcamp/seqprg/hisat2-2.0.4/hisat2 -p 4 --no-unal --dta -x /data/bootcamp/refdb/genome -1 SRR1551069_1.fastq.gz -2 SRR1551069_2.fastq.gz -S SRR1551069.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools view -b -u -S -o SRR1551069.unsort.bam SRR1551069.sam
/data/bootcamp/seqprg/samtools-1.3.1/samtools sort -o SRR1551069.bam SRR1551069.unsort.bam
/data/bootcamp/seqprg/samtools-1.3.1/samtools index SRR1551069.bam
`

While the alignment is running, we will examine the DNA alignment

If you are using windows, use WINSCP or your favorite data transfer tool to download the fastqc output onto your computer.
If you are using a mac, go to the terminal window where you previously downloaded files to your computer

`
scp username@toxea.swmed.edu:~/day1/NA12878.ba\* .
`

You can open the BAM file in IGV.  Open IGV -- Load the Genome using: Genomes->Load from Server->Select Human (b37)
Next upload the BAM files using File->Load From File

Try to answer the following questions:

* What is the range of coveage of the exons in CRLF2 *hint put your mouse on the histogram plot
* Look at the region: X:39,933,319-39,933,359

    * Can you see an SNV?
    * At what position?
    * How many reads supporting the alt allele?
    * How many reads supporting the ref allele?
    * What is the gene
    * What is the amino acid change? *Hint use google

* Look at the region: X:15,818,094-15,818,134

    * Can you see an indel?
    * At what positions?
    * How many reads supporting the alt allele?
    * How many reads supporting the ref allele?
    * What is the gene?
    * What is the indel an insertion or deletion?
    * What is it in an exon or intron?

* If you have time, examine these other positions:

    * X:15838366
    * X:41094888
    * X:41208310
    * X:70338399
    * X:70341169
    * X:70347330
    * X:70349947
    * X:100604757
    * X:123480147
    * X:129147079
    * X:129190217
    * X:152807923
    * X:152815089
    * X:152821717
    * X:152821887
    * X:152823728
    * X:152825414
    * X:152825433
    * X:152847291
    * X:152847301
    * X:152847488
    * X:152847553
    * X:152916023
    * X:152916236
    * X:153630255
    * X:153630293

Now it's time to examine the RNA alignments

If you are using windows, use WINSCP or your favorite data transfer tool to download the fastqc output onto your computer.
If you are using a mac, go to the terminal window where you previously downloaded files to your computer

`
scp username@toxea.swmed.edu:~/day1/SRR\*.ba\* .
`

You can open the BAM file in IGV.  Open IGV -- Load the Genome using: Genomes->Load from Server->Select Human (hg38)
Next upload the BAM files using File->Load From File

Please examine the following genes and type to answer these questions:
ALPL,CCNJL,KAZN,OLFM4,FCGR3B,PRSS21,PRMT6,SPOCK1,NEXN,ME1

* Do the reads cover all of the exons in the gene?
* Which sample has higher abundance?
* Are the exons covered equally?
* Are there differences in splice junctions
  * Hint Try Sashimi Plots:https://www.broadinstitute.org/igv/Sashimi

