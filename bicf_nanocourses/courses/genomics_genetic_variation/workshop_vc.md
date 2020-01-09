# Workshop for Variant Calling

Today we are going to:

* Run a germline variant calling program called strelka2 on 2 samples
* Run a somatic variant calling program called shimmer on 2 samples
* Run a structural variant calling program called delly on 2 samples
* Examine Alignments and variants in IGV and iobio
* Filter Variant for Quality

## Log into BioHPC
First, we will log into the log into a compute node,install a necessary package and then log into a compute node

* Set up a [WebGui](https://portal.biohpc.swmed.edu/terminal/webgui) session on BioHPC
* Log into via the VNC
* Open a terminal window -- you should be in the directory /archive/nanocourse/genome_analysis/trainXX (train01 used in the example)
* Copy data from /archive/nanocourse/genome_analysis/shared/vc_calling_session2 to your train directory
* hint 

~~~~
cd /archive/nanocourse/genome_analysis/shared/vc_calling_session2
cp * /archive/nanocourse/genome_analysis/train01/
cd /archive/nanocourse/genome_analysis/train01/vc_calling_session2
~~~~

#### Germline SNV calling

## Run Germline Variant Calling Program Strelka 2

First we are going to run variant calling using a program called Stelka2 on 2 engineered cell line mixtures available from Horizon Genomics, [HD728](https://www.horizondiscovery.com/tru-q-1-5-tier-reference-standard-hd728) and [HD752](https://www.horizondiscovery.com/tru-q-0-100-wildtype-reference-standard-hd752).  These sample were sequenced to ~1000X coverage.

Here are the commands that you will need to run Strelka2

~~~~
module load strelka/2.9.0 samtools/1.6
mkdir manta strelka
configureStrelkaGermlineWorkflow.py --bam HD752.nanocourse.bam --bam HD728.nanocourse.bam --referenceFasta /project/shared/bicf_workflow_ref/GRCh38/genome.fa --runDir strelka
strelka/runWorkflow.py -m local -j 8
~~~~


## Run Somatic Variant Calling Program Shimmer

#### Somatic variants

~~~~
module load shimmer/0.1.1
mkdir shimmer
shimmer.pl --minqual 25 --ref /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD752.nanocourse.bam HD728.nanocourse.bam --outdir shimmer 2> shimmer.err
cd shimmer
module load snpeff/4.3q
java -jar $SNPEFF_HOME/snpEff.jar GRCh38.86 somatic_diffs.vcf > somatic_diffs_annotate.vcf
~~~~


## Run Structural Variant Calling Program Delly

~~~~
module load delly2/v0.7.7-multi
delly2 call -t BND -o delly_translocations.bcf -q 30 -g /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD728.nanocourse.bam HD752.nanocourse.bam
delly2 call -t DUP -o delly_duplications.bcf -q 30 -g /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD728.nanocourse.bam HD752.nanocourse.bam
delly2 call -t INV -o delly_inversions.bcf -q 30 -g /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD728.nanocourse.bam HD752.nanocourse.bam
delly2 call -t DEL -o delly_deletion.bcf -q 30 -g /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD728.nanocourse.bam HD752.nanocourse.bam
delly2 call -t INS -o delly_insertion.bcf -q 30 -g /project/shared/bicf_workflow_ref/GRCh38/genome.fa HD728.nanocourse.bam HD752.nanocourse.bam
~~~~

samples.tsv - change normal to control 

~~~~
delly2 filter -t BND -o  delly_tra.bcf -f somatic -s samples.tsv delly_translocations.bcf`
delly2 filter -t DUP -o  delly_dup.bcf -f somatic -s samples.tsv delly_duplications.bcf`
delly2 filter -t INV -o  delly_inv.bcf -f somatic -s samples.tsv delly_inversions.bcf`
delly2 filter -t DEL -o  delly_del.bcf -f somatic -s samples.tsv delly_deletion.bcf`
delly2 filter -t INS -o  delly_ins.bcf -f somatic -s samples.tsv delly_insertion.bcf`

module load vcftools/0.1.14 samtools/1.6
bcftools concat -a -O v delly_dup.bcf delly_inv.bcf delly_tra.bcf delly_del.bcf delly_ins.bcf| vcf-sort -t temp > delly.vcf
java -jar $SNPEFF_HOME/snpEff.jar GRCh37.75 delly.vcf > delly_annotate.vcf
mkdir delly
mv delly_* delly
~~~~

## Run the Germline Program on Astrocyte

While you are running the pipeline on the command line, you can run the full analysis pipelines using our point and click workflows on bioHPC!

1. Go to [Astrocyte](https://astrocyte.biohpc.swmed.edu)
2. Create a new project test
3. Add data to your project including fastq and design file
4. Start Variant Germline Workflow using the Fastq files and the germline.design.txt file



## Examine the HD728.nanocourse.bam File in IGV

~~~~
module load IGV/2.3.90
sh igv.sh 
~~~~

You can open the BAM file in IGV.  Open IGV -- Load the Genome using: Genomes->Load from Server->Select Human (hg38)
Next upload the BAM files using File->Load From File

Try to answer the following questions:

* What is the range of coveage of the exons in KRAS *hint put your mouse on the histogram plot
* Look at the positions and answer the questions below:
    * chr7:140753336
    * chr7:140753337
    * chr7:55174014
    * chr7:55181378
    * chr13:28018500
    * chr2:208248389
    * chr9:5073770
    * chr12:25245350
    * chr12:25245351
    * chr12:25245347
    * chr15:66436825
    * chr9:136504914
    * chr1:114713909
    * chr3:179234297

* Can you see an SNV or Indel in these positions?
* At which position?
* How many reads supporting the alt allele?
* How many reads supporting the ref allele?
* What is the gene?
* What is the amino acid change? *Hint use google


## Load VCF into IGV

Compare your visual results to the variants identified manually from these position by loading the VCF file into IGV

* Which variants were called by the variant calling programs?
* Which variants are somatic and which are germline?

## Load BAMs and Variants into Iobio

* Load BAM to get [BAM stats](bam.iobio.io)
  * Select BAM and BAI file
* Load VCF to get [VCF stats](vcf.iobio.io)
  * Select VCF and VCF.TBI file (strelka2)
* Load BAM and VCF to get [Variant stats](gene.iobio.io)
  * Select BAM and BAI file
  * Select VCF and VCF.TBI file (strelka2)
  * Input the Genes that you found in the section above (hint click genes) and click analyze all
  * Pick a few genes (bookmark)and save to a file
  * Try filters and set minimum coverage > 25, allele frequency < 0.05 and exclude non-coding intronic variants


