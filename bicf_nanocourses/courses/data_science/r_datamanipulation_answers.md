# R Data Manipulation

## Height and Weight

In this exercise, we are going to explore a dataset of height and weight of teens

To load this data:
```
tbl <- read.csv(file='heightweight.csv')
```

To explore this data, you can either use the dplyr package (which has a variety of useful methods to filter, aggregate, and sort data frames), or you can try the sqldf package, which lets you run SQL queries on top of data frames, for even greater flexibility -- if SQL is your thing.

Explore the data using the tools you have. You can use the summary function, plot some of the data or look at histograms, and so on.  Try the following questions:

* Calculate the weight in Kg
```
tbl <- mutate(tbl,weightKg=weightLb/2.2)
```
* Calculate the  height in meters
```
tbl <- mutate(tbl,heightm=heightIn*0.0254)
```
* Calculate BMI as kg/m2
```
tbl <- mutate(tbl,bmi=weightKg/(heightm*heightm))
```
* Plot the BMI by gender using a box plot
```
boxplot(tbl$bmi ~ tbl$sex)
```
* What is the min, max, mean and median height per gender?
```
> summarize(group_by(tbl,sex),mean.height=mean(heightIn),min.height=min(heightIn),max.height=max(heightIn),median.height=median(heightIn))
# A tibble: 2 x 5
     sex mean.height min.height max.height median.height
  <fctr>       <dbl>      <dbl>      <dbl>         <dbl>
1      f    60.52613       51.3       66.8          61.3
2      m    62.10317       50.5       72.0          61.9
```

## Bacterial Community

In this exercise, we are going to explore the bacterial communities of 4 synthetic samples.  This dataset contains the counts of reads that are classified into each taxa per sample.  The first column is the taxonomic level where:

1. Kingdom
2. Phylum
3. Class
4. Order
5. Family
6. Genus
7. Species

To load the data use the command:
```
mock <-read.table(file="mock_community.txt",sep='\t',header=TRUE)
```

To explore this data, you can either use the dplyr package (which has a variety of useful methods to filter, aggregate, and sort data frames), or you can try the sqldf package, which lets you run SQL queries on top of data frames, for even greater flexibility -- if SQL is your thing.

Explore the data using the tools you have. You can use the summary function, plot some of the data or look at histograms, and so on.  Try the following questions:

* What is the most abundant genus for each community?
```
> sqldf("select taxon,Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=6 order by Mock1 desc limit 1")
           taxon Mock1 Mock2 Mock3 Mock4
1 Staphylococcus    68  1327   256   295
> sqldf("select taxon,Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=6 order by Mock2 desc limit 1")
                    taxon Mock1 Mock2 Mock3 Mock4
1 Bacillales_unclassified    50  1472   380   446
> sqldf("select taxon,Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=6 order by Mock3 desc limit 1")
                    taxon Mock1 Mock2 Mock3 Mock4
1 Bacillales_unclassified    50  1472   380   446
> sqldf("select taxon,Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=6 order by Mock4 desc limit 1")
                    taxon Mock1 Mock2 Mock3 Mock4
1 Bacillales_unclassified    50  1472   380   446
```
* Which phylum has the highest read counts for all 4 samples?  Can you create an sql statement to determine the phylum with the max sum?
```
> sqldf("select taxon,Mock1+Mock2+Mock3+Mock4 as TotalAbund from mock where taxlevel=2 order by TotalAbund desc limit 1")
       taxon TotalAbund
1 Firmicutes       7897
```
* Calculate the averge number of reads for each phyla
```
 sqldf("select taxon,(Mock1+Mock2+Mock3+Mock4)/4 as TotalAbund from mock where taxlevel=2")
```

* Create a new table with the sample name and the ratio of bacteroidetes/firmicutes.
```
taxnames <- as.character(sqldf("select taxon from mock where taxlevel=2")[,1])
phylum <- t(sqldf("select Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=2"))
colnames(phylum) <-taxnames
mutate(as.data.frame(phylum),ratio=Bacteroidetes/Firmicutes)
```
* Which sample has the highest ratio?
```
Mock1
```
* Create a table with just the genera (genus ie taxlevel=6).  Add the total number of reads per genera and sort the table by the total.
```
> genera <- sqldf("select taxon,Mock1,Mock2,Mock3,Mock4,Mock1+Mock2+Mock3+Mock4 as TotalAbund from mock where taxlevel=6 order by TotalAbund")
```
* Create a heatmap or stack bar plot of the genera by sample.
```
library(reshape)
library(ggplot)
phylum <- sqldf("select taxon,Mock1,Mock2,Mock3,Mock4 from mock where taxlevel=6")
taxct <- melt(phylum,id.vars=c('taxon'))
names(taxct) <- c('taxname','sample','taxct')
total <- sqldf("select sample,sum(taxct) as nseqs from taxct group by sample")
taxct <- merge(taxct,total,by='sample')
taxct <- mutate(taxct,perc=round(taxct/nseqs, digits = 7))
min.seq<- min(total$nseqs)
taxct <- mutate(taxct,norm=min.seq*perc)
ggplot(taxct, aes(x = sample, y = perc,fill=taxname)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Now load a second table of "real" samples:
```
real <-read.table(file="real_community.txt",sep='\t',header=TRUE)
```

* Join mock with real
```
join.tbl <- sqldf("select mock.taxlevel,mock.taxon,Mock1,Mock2,Mock3,Mock4,Sample1,Sample2 from mock inner join real using(taxon)")
```

* Create a table at any taxa level (ie phylum, order, etc) and plot taxonomic abundance at of samples with stack barplots
```
library(reshape)
library(ggplot)
taxct <- melt(join.tbl,id.vars=c('taxon','taxlevel'))
names(taxct) <- c('taxname','level','sample','taxct')
total <- sqldf("select sample,sum(taxct) as nseqs from taxct group by sample")
taxct <- merge(tax,total,by='sample')
taxct <- mutate(taxct,perc=round(taxct/nseqs, digits = 7))
min.seq<- min(total$nseqs)
taxct <- mutate(taxct,norm=min.seq*perc)
ggplot(taxct, aes(x = sample, y = perc,fill=taxname)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
```




