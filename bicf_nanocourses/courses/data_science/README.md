# [BICF](http://www.utsouthwestern.edu/labs/bioinformatics/) Data Science for Biologist

Do you want to be to do simple statistical analyses yourself? Do you find yourself spending time and effort generating the same plots and statistics for each project? R is a freely available language and programming environment for statistical computing and graphics which provides a wide variety of statistical and graphical techniques: linear and nonlinear modelling, statistical tests, time series analysis, classification, clustering, etc.

[Download R Studio](https://www.rstudio.com/products/rstudio/download/)

If you already have R studio, please [update your R here](https://cran.r-project.org/)

**Here are datasets that you will need for this course:**

[Statistical Test Dataset](dig_csv.zip)<br>
[RNA-seq Dataset](countTable.txt)

***
## Contacts
* Course Coordinator [Brandi Cantarel](mailto:brandi.cantarel@utsouthwestern.edu)
* Course Administration [Neha Sinha](Neha.Sinha@UTSouthwestern.edu)
* TAs
    - Mingzhu Nie (Jan 10)
    - Erika Villa (Jan 10)
    - Holly Ruess (Jan 17)
    - Krishna Kanth Chitta (Jan 17)
    - Jaideep Chaudhary (Jan 24 and Jan 31)

***

## R Cheatsheets

* [Data Import](data-import.pdf)
* [String Functions](strings.pdf)
* [Dates and Times](lubridate.pdf)
* [More Cheatsheets](https://rstudio.com/resources/cheatsheets/)

## Schedule

|  Topic | Instructor|
| ------------- | ------------- |
| 1/10/2020 | Room NB2.100A|
| [Introduction to R](RDataStructures.pdf) and [R Data Structures](r_intro.html)<br>[Rmd](r_intro.Rmd)<br>[Workshop 1](rdatastructures.md)<br>[DataFiles](DataStructureLecture.zip) | Brandi Cantarel |
| [Data Importing and Cleaning with Tidyverse](DataClean.html)<br>[Workshop 2](DataCleaning.zip) | Brandi Cantarel |
| [Data Manipulation and Data Joining with dplyr and tidyverse](r_dataManipulation_lecture_1_8_20.zip)<br>[Workshop]( r_dataManipulation_workshop_1_8_20.zip)<br>[Workshop Solutions](dataManipulation_workshop_answers.html)| Spencer Barnes |
| 1/17/2020 | Room NB2.100A|
| [Introduction to Statistical Tests](Nanocourse_StatisticalTests.html)<br>[Workshop](Nanocourse_StatisticalTests_Workshop_Blank.Rmd)<br>[Workshop Solutions](Nanocourse_StatisticalTests_Workshop.html) | Jeremy Mathews |
| [Correlations and Linear Regression](Nanocourse_CorrelationLR.html)<br>[Workshop](Nanocourse_CorrelationLR_Workshop_Blank.Rmd)<br>[Workshop Solutions](Nanocourse_CorrelationLR_Workshop.html)|  Jeremy Mathews |
| [Plotting with GGPlot and plotly](Plotting_with_GGplot_and_Plotly.nb.html) <br>[Data](Data_PlottingWithGGplot_Ploty.zip) <br>[Workshop](Workshop_Plotting_with_GGplot_and_Plotly.nb.html) | Jeon Lee |
| 1/24/2020 | Room NB2.100A|
| [Programming basics](programing_basics_cheatsheet.md) <br> [Psudeo Code](programing_basics_psudocode.md) | Venkat Malladi |
| [Loops and Looping functions with Apply, Scripting and Markdown](ScriptingMarkdown_lecture.html)<br>[R code](ScriptingMarkdown_lecture.Rmd)<br>[Workshop](rscripting.Rmd) <br>[Data](Nanocourse_R_II_Jan2020.zip) | Chris Bennett |
| 1/31/2020 | Room NG3.202 |
| R Package Repositories<br>[Presentation](PackageRepositories.pptx)<br>[Workshop](PackageRepositories.R)<br>[RNA-seq Dataset](countTable.txt)<br>[Workshop Answers](PackageRepositories_answers.R) | Gervaise Henry |
| Accessing Public Data Though Bioconductor<br>[TCGAbiolinks Lecture](tcga_nanocourse_lecture.nb.html)<br>[TCGAbiolinks Workshop](tcga_workshop_questions.nb.html)<br>[Workshop Answers](tcga_workshop_answers.nb.html) | Spencer Barnes |
| Student Projects | |

***
### Student Projects
Here the an opportunity to apply what you learned to your own research!  Students should present a question with some possible solutions to discuss as a class.

Here are some example questions:

1. Pick a dataset from the class or from your own work.  Plot 2 continuous variables and add a trend lines.  Create box-lots using a continuous variable and a categorical variable.  Add text to indicate the mean of each group on the plot (type mtext)  Present the summary statistics for the comparison between a continuous variable and a categorical variable.
2. Calculate logCMP from RNASeq read count data and make a heatmap of the a subset of genes — chose 10 or 20 genes.  Choose 2 genes to make boxplots comparing the expression with the sample groups.  Create a 3-D plots to show the expression of 3 genes.
3. Create 3 vectors using random functions for classic distributions see  https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html Plot these vectors as histograms, cumulative distribution function and density function.
4. Pick a package in bioconductor, prepare 5-10 slides to show the other students in the class on how to use this function from installation to some final plot.  
5. Pick a plot from a recent publication and determine how to make that plot in R.
