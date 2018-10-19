
# Workshop for Variant Prioritization

Today we are going to:
 
* Practice VCF manipulation skills with BCFtools
* Perform project QC
* Filter variants with GEMINI
    * Solving a rare disease case from a trio.
    * Identifying somatic muations in paired cancer-normal cell lines.

## Log into BioHPC
First, we will log into the a compute node

* Set up a [WebGUI](https://portal.biohpc.swmed.edu/terminal/webgui) session on BioHPC
* Launch via "connect with VNC client", open using [turbovnc](https://sourceforge.net/projects/turbovnc/). You can also launch the session by "connect via web" but copying and pasting may not work under this mode.
* Open a terminal window -- you should be in the directory /archive/nanocourse/May2018/trainXX
* Copy session3 material into your directory and work from there
~~~~
cp -r /archive/nanocourse/May2018/shared/session3 .
cd session3
~~~~


## Practice VCF manipulation skills with BCFtools
First we are going to practice VCF manipulation skills on VCF from [CEPH family 1463 with 17 members](https://www.coriell.org/0/Sections/Collections/NIGMS/CEPHFamiliesDetail.aspx?PgId=441&fam=1463&).

* Load bcftools module on BioHPC
~~~~
module load bcftools 
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
tabix -l vcf_playground/maternal_family_recombined.vcf.gz
~~~~



## Perform project QC

* Run pedigree VCF QC tool [peddy](https://github.com/brentp/peddy)
~~~~
source activate py2.7-nano
mkdir peddy_QC
cd peddy_QC
python -m peddy -p 16 --plot --prefix ceph-1463 ../ceph1463.vcf.gz ../ceph1463.ped
cd ..
~~~~

## Using GEMINI

~~~~
module load gemini
mkdir gemini_practices
mkdir ~/.gemini
cp gemini-config.yaml ~/.gemini/
~~~~

* Look at gemini usage messages  
~~~~
gemini --help
~~~~
~~~~
gemini load --help
gemini annotate --help
gemini query --help
~~~~
~~~~
gemini de_novo --help
gemini comp_hets --help
gemini autosomal_recessive --help
gemini x_linked_recessive --help
~~~~
~~~~
gemini set_somatic --help
~~~~
We will try out some of these tools in the following commands, you may refer to the documentation to understand the options we will be using.

### Solving a rare disease case from a trio

#### A patient presented at the hospital with hyperammonemia after giving birth. The patient and her parents have been sequenced. The VCF (Hg19) and ped files are provided for you.  Can you find out the disease causing variant(s)?
* Load vcf and ped files into a GEMINI database **(this can take a while, you may skip this step and proceed from GMDP_3_0054.db, which is the result from this step that has already been generated for you)**.
~~~~
gemini load -v GMDP_3_0054.annot.vt.vep.vcf.gz -t VEP -p GMDP_3_0054.ped --cores 32 GMDP_3_0054.db
~~~~
* Since the vcf was generated as a union of four callers, we want to import the information from "CallSet" field of the vcf into the gemini database as well.
~~~~
gemini annotate -f GMDP_3_0054.annot.vt.vep.vcf.gz \
    -a extract -c CallSet -t text -e CallSet -o uniq_list GMDP_3_0054.db
~~~~
* What are the tables and fields for each in this database?
~~~~
gemini db_info GMDP_3_0054.db
~~~~
* Take a look at the ped file
~~~~
cat GMDP_3_0054.ped
~~~~
* Check the sample names
~~~~
gemini query -q "SELECT name FROM samples" --header GMDP_3_0054.db
~~~~
* Which sample is from the affected proband (phenotype being 2)?
~~~~
gemini query -q "SELECT name FROM samples WHERE phenotype == 2" GMDP_3_0054.db
~~~~
* How many variants are from chromosome 1?
~~~~
gemini query -q "SELECT COUNT(*) FROM variants WHERE chrom == 'chr1'" GMDP_3_0054.db
~~~~

#### Pathogenic variants from clinvar
~~~~
gemini query -q \
    "select gene,chrom,start,end,ref,alt,impact,codon_change,aa_change,\
    gt_types,gt_depths,gt_ref_depths,gt_alt_depths,max_aaf_all,gnomad_num_hom_alt,gnomad_num_het,\
    clinvar_sig,clinvar_disease_name,clinvar_gene_phenotype,gerp_bp_score,cadd_scaled \
    from variants where clinvar_sig LIKE '%pathogenic%'" \
    --header GMDP_3_0054.db > gemini_practices/GMDP_3_0054.pathogenic.txt
~~~~

#### Filter by different inheritance models
* de novo 
~~~~
gemini de_novo GMDP_3_0054.db \
    --columns "gene,chrom,start,end,ref,alt,impact,codon_change,aa_change,\
    gt_types,gt_depths,gt_ref_depths,gt_alt_depths,max_aaf_all,gnomad_num_hom_alt,gnomad_num_het,\
    clinvar_sig,clinvar_disease_name,clinvar_gene_phenotype,gerp_bp_score,cadd_scaled" \
    --filter "max_aaf_all < 0.001 AND impact_severity != 'LOW' AND gnomad_num_het < 1">\
    gemini_practices/GMDP_3_0054.denovo.txt
~~~~
* compound heterozygous
~~~~
gemini comp_hets GMDP_3_0054.db \
    --columns "gene,chrom,start,end,ref,alt,impact,codon_change,aa_change,\
    gt_types,gt_depths,gt_ref_depths,gt_alt_depths,max_aaf_all,gnomad_num_hom_alt,gnomad_num_het,\
    clinvar_sig,clinvar_disease_name,clinvar_gene_phenotype,gerp_bp_score,cadd_scaled" \
    --filter "max_aaf_all < 0.001 AND impact_severity != 'LOW' AND gnomad_num_het < 1">\
    gemini_practices/GMDP_3_0054.comp_hets.txt
~~~~
* autosomal recessive
~~~~
gemini autosomal_recessive GMDP_3_0054.db \
    --columns "gene,chrom,start,end,ref,alt,impact,codon_change,aa_change,\
    gt_types,gt_depths,gt_ref_depths,gt_alt_depths,max_aaf_all,gnomad_num_hom_alt,gnomad_num_het,\
    clinvar_sig,clinvar_disease_name,clinvar_gene_phenotype,gerp_bp_score,cadd_scaled" \
    --filter "max_aaf_all < 0.001 AND impact_severity != 'LOW' AND gnomad_num_het < 1">\
    gemini_practices/GMDP_3_0054.autosomal_recessive.txt
~~~~
* x-linked recessive
~~~~
gemini x_linked_recessive GMDP_3_0054.db \
    --columns "gene,chrom,start,end,ref,alt,impact,codon_change,aa_change,\
    gt_types,gt_depths,gt_ref_depths,gt_alt_depths,max_aaf_all,gnomad_num_hom_alt,gnomad_num_het,\
    clinvar_sig,clinvar_disease_name,clinvar_gene_phenotype,gerp_bp_score,cadd_scaled" \
    --filter "max_aaf_all < 0.001 AND impact_severity != 'LOW' AND gnomad_num_het < 1">\
    gemini_practices/GMDP_3_0054.x_linked_recessive.txt
~~~~

#### Number of records under each inheritance mode
~~~~
for file in gemini_practices/GMDP_3_0054.*.txt; do wc -l $file; done
~~~~
You may open these files in excel and check further

#### Making a diagnosis
* Find out genes and HPO terms associated with Hyperammonemia using resource downloaded from [here](http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/)
* You may also search for the term [here](http://compbio.charite.de/hpoweb/)
~~~~
grep Hyperammonemia ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt > gemini_practices/Hyperammonemia_candidate_genes.txt
~~~~
* Alternatively you can use [Phenomizer](http://compbio.charite.de/phenomizer/) to obtain a list of candidate diseases and related gene sets. This would be particularly helpful when multiple phenotypes present.

Can you find out what might be the causal variant?  

* Answer
~~~~
cd gemini_practices
# for each file generated by a specific inheritance model or ad hoc query:
for file in GMDP_3_0054.*.txt; do 
    # select unique genes from the first column as our candidate genes, loop through them:
    for candidate_gene in $(cut -f1 $file|sort|uniq); do 
        # for each known disease gene in the second column of the file we obtained through HPO term mapping:
        for disease_gene in $(cut -f2 Hyperammonemia_candidate_genes.txt); do 
            # match to see if the candidate gene is in a known disease gene, if it is, then
            if [ $candidate_gene = $disease_gene ]; then 
                # print the name of the disease gene, and the model name of the file
                echo "$disease_gene found under $(echo $file|cut -f2 -d .) mode"
                # create a result file (if it does not exist)
                touch trio_result.txt
                # write the header into this result file
                head -n1 $file >> trio_result.txt
                # append the variant records to this result file
                grep $disease_gene $file >> trio_result.txt
            fi
        done
    done
done
~~~~
~~~~
# check columns 1-6 and 10 of this result file
cut -f -6,10 trio_result.txt | column -t | less -S
~~~~

### Identifying somatic muations in paired cancer-normal cell lines.

#### A pair of cell lines from the same lung cancer patient (HCC4017, cancer cell line vs HBEC30, immortalized human bronchial epithelial cell line) have been sequenced and the VCF (lifted over from Hg38 to Hg19) and ped files are provided to you. Can you identify the somatic mutations?
* Load vcf and ped files into a GEMINI database **(this can take a while, you may skip this step and proceed from CellLine.db, which is the result from this step that has already been generated for you)**.
~~~~
gemini load -v CellLine.vt.vcf.gz -t snpEff -p CellLine.ped --cores 32 CellLine.db
~~~~
* Check the ped file
~~~~
cat CellLine.ped
~~~~
* Check the sample names
~~~~
gemini query -q "SELECT name FROM samples" --header CellLine.db
~~~~
* Find run set_somatic to flag somatic mutations and save output
~~~~
(cat somatic_mutation_header; \
gemini set_somatic \
    --min-depth 30 \
    --min-qual 20 \
    --min-somatic-score 15 \
    --min-tumor-depth 10 \
    --min-norm-depth 10 \
    CellLine.db) > gemini_practices/CellLine_Somatic_Mutations.txt
~~~~

* Obtain additional information for somatic mutations  
~~~~
gemini query -q \
    "SELECT chrom,start,end,gene,vcf_id,somatic_score,gt_types,ref,alt,qual,type,sub_type,codon_change,aa_change,\
    aa_length,biotype,impact_severity,impact,is_somatic,gerp_bp_score,cadd_scaled,fitcons\
    FROM variants WHERE is_somatic==1 ORDER BY vcf_id DESC" \
    --header CellLine.db > gemini_practices/CellLine_Somatic_Mutations_additional.txt
~~~~
~~~~
column -t gemini_practices/CellLine_Somatic_Mutations_additional.txt | less -S
~~~~
You may open these files in excel and check further
