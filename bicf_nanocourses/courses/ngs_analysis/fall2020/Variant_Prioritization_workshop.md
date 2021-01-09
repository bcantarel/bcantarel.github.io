# Workshop for Variant Prioritization

Today we are going to:
 

* Find genes that are associated with a phenotype
* Identify Candidate variation
    * Solving a rare disease case from a trio.
    * Identifying somatic muations in paired cancer-normal cell lines.
* Filter Varinats by QC using BCFTools

## Data for this Workshop

* [To mount bioHPC storage onto your computer](https://portal.biohpc.swmed.edu/content/guides/biohpc-cloud-storage/)
* [Lamella Data Storage](https://lamella.biohpc.swmed.edu/index.php/apps/files/)
  * Click your profile icon on the right and select settings
  * Select External Storage on the Left
  * Add a new storage called course (on left)
  * Select "Username and Password"
  * Add on the right: "archive", "nanocourse/sequence_analysis", "trainXX", "trainX_password"
  
* Navigate to /archive/nanocourse/sequence_analysis/shared/variant_prioritization_session2


## Find Genes associated with a Phenoype using the NCBI Tools

### Find Variants of a gene of Interest
* From the [NCBI Gene database](https://www.ncbi.nlm.nih.gov/gene)
1. Search the Term = msh6 (also try this in PubMed) or msh6[sym] AND human[orgn]
2. Go to the Variation section of the Gene record (in the Table of Contents)
3. Variation Viewer (GRCh38) 
4. Filter by Source dbSNP, in ClinVar, Pathogenic

### Find Variants Associated with a Disease
* From the [MedGen database](https://www.ncbi.nlm.nih.gov/medgen)
1. Search the term = severe combined immunodeficiency
2. Some have active Gene links. To limit to those that do:
 - Use Limits (link under search box) to select those "Associated with a gene", run Search
 - Use the "Find related data" menu, select the Gene database (takes you to Gene)
3. Pick of of the listed diseases and click the ClinVar link (open in new window)
 - Filter Pathogenic Variants
 - Filter by Reviewer Status to find expertly reviewed variants.
 - How many variants do you find
4. Select on link for the disease
    - What are tissues where the genes associated with that disease are expressed?

## Filtering Variants

### Solving a rare disease case from a trio

#### Data 
Find data files here: /archive/nanocourse/sequence_analysis/shared/variant_prioritization_session2

#### A patient presented at the hospital with hyperammonemia after giving birth. The patient and her parents have been sequenced. The VCF (Hg19) file provided for you.

* Go to [Gene.iobio](http://gene.iobio.io)
* Load the VCF File and Index, hyperammonemia.vcf.gz and hyperammonemia.vcf.gz.tbi by Clicking Files
   - Proband = GMDP_3_0054_1
   - Father = GMDP_3_0054_2
   - Mother = GMDP_3_0054_3
   - Hint click trio and load the file 3 files and select the sample name for each of family member
* Add a gene list
   - Use the phenotype button to add the gene list.
   - Increase the number of genes from 10 to 100
* See all variants (not just coding variants) -- Under Options
* Also update filters to see all impact variants not just those with "know pathgenicity"
* Once you have updated the filters, click Analyze all

Hint, this is a compound het with one mutation not in the coding region

Can you find out what might be the causal variant?  

### Identifying disease causing muations in paired cancer sample

#### Data
 Find data files here: /archive/nanocourse/sequence_analysis/shared/variant_prioritization_session2

#### A patient presented at the hospital with .  The patient's tumor has been sequenced.  The VCF (HG38) file has been provided for you.

* Here is a gene list for clinically actionable solid tumors

* Go to [Gene.iobio](http://gene.iobio.io)
* Load the VCF File and Index, AML_Cancer.vcf.gz and AML_Cancer.vcf.gz.tbi
* Change the Genome Build to HG38
* Add a gene list
   - AKT1,ALK,APC,ATM,BRAF,CDH1,CTNNB1,EGFR,ERBB2,FBXW7,FGFR2,GNAQ,GNAS,HRAS,IDH1,KIT,KRAS,MET,NRAS,PDGFRA,PIK3CA,PTEN,SMAD4,SMARCB1,SMO,SRC,STK11,TP53
* Explore these variants
   - Are any of these Pathogenic or High Impact variants seen in other AML patients (search COSMIC)
   - Are any of these mutations associated with treatment (search CIVIC)
   - Are any of these mutations seen in subjects in GnomAD
 
## Log into BioHPC
First, we will log into the a compute node

* Set up a [WebGUI](https://portal.biohpc.swmed.edu/terminal/webgui) session on BioHPC
* Launch via "connect with VNC client", open using [turbovnc](https://sourceforge.net/projects/turbovnc/). You can also launch the session by "connect via web" but copying and pasting may not work under this mode.
* Open a terminal window -- you should be in the directory /archive/nanocourse/genome_analysis/trainXX
* Copy session3 material into your directory and work from there
~~~~
cp -r /archive/nanocourse/sequence_analysis/shared/variant_prioritization_session2 .
cd variant_prioritization_session2
~~~~


## Practice VCF manipulation skills with BCFtools
We are going to practice VCF manipulation skills on VCF from [CEPH family 1463 with 17 members](https://www.coriell.org/0/Sections/Collections/NIGMS/CEPHFamiliesDetail.aspx?PgId=441&fam=1463&).

* Load bcftools module on BioHPC
~~~~
module load bcftools htslib/gcc/1.8
~~~~


* Compress the VCF
~~~~
bgzip -c ceph1463.vcf > ceph1463.vcf.gz
~~~~
The "-c" option write on standard output, keep original files unchanged.  
If you don't want to keep the original file, do:
~~~~
bgzip ceph1463.vcf
~~~~
If you want to decompress do
~~~~
bgzip -d ceph1463.vcf.gz
~~~~

* Generate index (.tbi file) using tabix (loaded with the bcftools module)
~~~~
tabix ceph1463.vcf.gz
~~~~

* Build a new directory to practice VCF manipulation skills
~~~~
mkdir vcf_playground
~~~~

* Look at bcftools usage messages
~~~~
bcftools --help
bcftools query --help
bcftools stats --help
bcftools filter --help
bcftools view --help
~~~~
We will try out some of these tools in the following commands, you may refer to the documentation to understand the options we will be using.


* What are the samples in this VCF?
~~~~
bcftools query -l ceph1463.vcf.gz
~~~~

* Calculate stats on VCF, how many SNPs, MNPs and indels?  
~~~~
bcftools stats ceph1463.vcf.gz > vcf_playground/ceph1463.stats.out
less vcf_playground/ceph1463.stats.out
~~~~

* Extract just the chromosome, position and genotypes
~~~~
bcftools query -f '%CHROM\t%POS\t[%GT ]\n' ceph1463.vcf.gz | less -S
~~~~

* Extract region with bcftools
~~~~
bcftools filter --targets 10:96447911-96613017 ceph1463.vcf.gz | less
~~~~

* How many variants in this VCF? (grep -v "^#" excludes the meta-info and header lines that start with "#")
~~~~
zcat ceph1463.vcf.gz | grep -v "^#" |wc -l
~~~~

* How many variants remaining after filtering out SNPs within 10bp of an indel?
~~~~
bcftools filter --SnpGap 10 ceph1463.vcf.gz | grep -v "^#" | wc -l
~~~~

* Create a subset VCF with the maternal family members (-O specifies the compression format, -o specifies output file)
~~~~
bcftools view --samples NA12891,NA12892,NA12878 ceph1463.vcf.gz -O z -o vcf_playground/maternal_family.vcf.gz
~~~~

* Check the samples in the new VCF
~~~~
bcftools query -l vcf_playground/maternal_family.vcf.gz
~~~~

* Create a subset VCF with just the indels (-i specifies the inclusion expression)
~~~~
bcftools filter -i 'TYPE="indel"' vcf_playground/maternal_family.vcf.gz -O z -o vcf_playground/maternal_indels.vcf.gz
~~~~

* Check if we now just have the indels in this subset VCF (grep -A n print n lines after matched pattern)
~~~~
bcftools stats vcf_playground/maternal_indels.vcf.gz | grep -A 8 "SN, Summary numbers:"
~~~~

* Some of the indels have genotype 0/0 in all members of the maternal family, suggesting they are from the paternal family.
~~~~
bcftools query -f '%CHROM\t%POS\t[%GT ]\n' vcf_playground/maternal_indels.vcf.gz | less -S
bcftools query -f '[%GT ]\n' vcf_playground/maternal_indels.vcf.gz|grep '0/0 0/0 0/0' | wc -l
~~~~

* Exclude such variants to keep only those from the maternal family.
~~~~
bcftools view --private -s NA12891,NA12892,NA12878 vcf_playground/maternal_indels.vcf.gz -O z -o vcf_playground/maternal_only_indels.vcf.gz
bcftools query -f '[%GT ]\n' vcf_playground/maternal_only_indels.vcf.gz | grep '0/0 0/0 0/0' | wc -l
~~~~

* Create a subset with the autosomes and a subset with the X chromosome and combine them back together
~~~~
tabix vcf_playground/maternal_family.vcf.gz
tabix -l vcf_playground/maternal_family.vcf.gz
~~~~
~~~~
bcftools view vcf_playground/maternal_family.vcf.gz --targets ^X -O z -o vcf_playground/maternal_family_autosomes.vcf.gz
tabix vcf_playground/maternal_family_autosomes.vcf.gz
tabix -l vcf_playground/maternal_family_autosomes.vcf.gz
~~~~
~~~~
bcftools view vcf_playground/maternal_family.vcf.gz --targets X -O z -o vcf_playground/maternal_family_X.vcf.gz
tabix vcf_playground/maternal_family_X.vcf.gz
tabix -l vcf_playground/maternal_family_X.vcf.gz
~~~~
~~~~
bcftools concat vcf_playground/maternal_family_autosomes.vcf.gz vcf_playground/maternal_family_X.vcf.gz -O z -o vcf_playground/maternal_family_recombined.vcf.gz
tabix vcf_playground/maternal_family_recombined.vcf.gz
tabix -l vcf_playground/maternal_family_recombined.vcf.gz
~~~~
