# R Data Manipulation

## Height and Weight

In this exercise, we are going to explore a dataset of height and weight of teens

To load this data:
```
tbl <- read.csv(file='heightweight.csv')
```

To explore this data, you can either use the dplyr package (which has a variety of useful methods to filter, aggregate, and sort data frames), or you can try the sqldf package, which lets you run SQL queries on top of data frames, for even greater flexibility -- if SQL is your thing.

Explore the data using the tools you have. You can use the summary function, plot some of the data or look at histograms, and so on.  Try the following questions:

1. Calculate the weight in Kg
2. Calculate the  height in meters
3. Calculate BMI as kg/m2
4. Plot the BMI by gender using a box plot
5. What is the min, max, mean and median height per gender?

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

1. What is the most abundant genus for each community?
2. Which phylum has the highest read counts for all 4 samples?  Can you create an sql statement to determine the phylum with the max sum?
4. Calculate the averge number of reads for each phyla
5. Create a new table with the sample name and the ratio of bacteroidetes/firmicutes.
6. Which sample has the highest ratio?
7. Create a table with just the genera (genus ie taxlevel=6).  Add the total number of reads per genera and sort the table by the total.
8. Create a heatmap or stack bar plot of the genera by sample.


Now load a second table of "real" samples:
```
real <-read.table(file="real_community.txt",sep='\t',header=TRUE)
```

9. Join mock with real
10. Create a table at any taxa level (ie phylum, order, etc) and plot taxonomic abundance at of samples with stack barplots




