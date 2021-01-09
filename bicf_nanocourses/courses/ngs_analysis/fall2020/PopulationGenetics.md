# Population Genetics Workshop

# Contents
* [Quality control of germline mutations](#quality-control)
* [Detect population stratification](#population-stratification)
* [Perform association test](#association-test)

# Prepare the environment and data
```
ln -s /archive/nanocourse/May2018/shared/train04/data/* .
export PATH=/archive/nanocourse/May2018/shared/train04/bin/:$PATH
module load R/3.4.1-gccmkl
alias ll='ls -lhtr'
```

# Quality control

1. Check the vcf file
```
zcat gws.vcf.gz | less -S
```

2. Format the file to plink format
```
plink --vcf gws.vcf.gz --recode --out gws
less -S gws.ped
less -S gws.map
```

3. Check minor allele frequency (MAF)
```
plink --file gws --freq --out gws
sort -k5,5g gws.frq | less
```

4. Keep only the SNPs with MAF > 0.01
```
plink --file gws --maf 0.01 --recode --out gws.maf_ge_0.01
plink --file gws.maf_0.01 --freq --out gws.maf_0.01
sort -k5,5g gws.maf_0.01.frq | less
```

5. Calculate the p-value of hardy weinberg equilibrium
```
plink --file gws.maf_0.01 --hardy --out gws.maf_0.01
sort -k9,9g gws.maf_0.01.hwe | less
```

6. Keep only the SNPs with p-value(HWE) > 0.001
```
plink --file gws.maf_0.01 --hwe 0.001 --recode --out gws.maf_0.01.hwe_0.001
```

# Population stratification

1. Combine GWS and 1000G
```
plink --file gws.maf_0.01.hwe_0.001 --merge 1000g --recode --out gws.maf_0.01.hwe_0.001.1000g
```

2. Calculate principal components of GWS and 1000 Genomes
```
plink --file gws.maf_0.01.hwe_0.001.1000g --pca --out gws.maf_0.01.hwe_0.001.1000g
less -S gws.maf_0.01.hwe_0.001.1000g.eigenvec
```

3. Combine the PCs and population groups
```
sort gws.maf_0.01.hwe_0.001.1000g.eigenvec > gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort
paste -d' ' gws.1000g.pop gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort | cut -f1,2,5- > gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort.pop
```

4. Plot the PCs for 1000G, GWS and 1000G+GWS
```
./plot_pc.sh gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort.pop gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort.pop.pdf
```

5. Exclude PC outliers
```
cat gws.maf_0.01.hwe_0.001.1000g.eigenvec.sort.pop | awk '{if($2=="GWS" && $3>0.01)print $1,$1}' > gws.maf_0.01.hwe_0.001.1000g.eigenvec.outlier
plink --file gws.maf_0.01.hwe_0.001 --remove gws.maf_0.01.hwe_0.001.1000g.eigenvec.outlier --recode --out gws.maf_0.01.hwe_0.001.no_pc_outlier
```

6. Create soft link to QC passed files
```
ln -s gws.maf_0.01.hwe_0.001.no_pc_outlier.ped gws.pass.ped
ln -s gws.maf_0.01.hwe_0.001.no_pc_outlier.map gws.pass.map
```

# Association test

1. Check phenotype file
```
head pheno.pass.ldl
head pheno.pass.age_sex
```

2. Re-generate PCs for QC passed data
```
plink --file gws.pass --pca --out gws.pass
```

3. Generate covariate file (age + sex + pc1 + pc2)
```
cat gws.pass.eigenvec | cut -f1-4 -d' ' > gws.pass.pc2
echo "FID IID pc1 pc2" > gws.pass.pc2.fmt
cat gws.pass.pc2 >> gws.pass.pc2.fmt
head gws.pass.pc2.fmt
paste -d' ' pheno.pass.age_sex gws.pass.pc2.fmt | cut -f1-4,7,8 -d' ' > pheno.pass.age_sex_pc2
head pheno.pass.age_sex_pc2
```

4. LDL association test (linear regression)
```
plink --file gws.pass --allow-no-sex --pheno pheno.pass.ldl --covar pheno.pass.age_sex --linear --out gws.pass.ldl
less gws.pass.ldl.assoc.linear
cat gws.pass.ldl.assoc.linear | awk '{if($5=="ADD" || $5=="TEST")print}' > gws.pass.ldl.assoc.linear.res
head gws.pass.ldl.assoc.linear.res
```

5. QQ plot
```
./plot_qq.sh gws.pass.ldl.assoc.linear.res gws.pass.ldl.assoc.linear.res.qqplot.pdf
```
