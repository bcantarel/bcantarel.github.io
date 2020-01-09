# Objects in R: Introduction to R Data Structures

## R Data Structures
* Vectors 
    - Atomic vector — a collection of values
    - Factors — special vectors that represent categorical data
* Matrix — a special vector with rows and columns
* Data frame — a special data structure of rows and columns, the default structure for reading in “excel-like” files
* List — a vector of different data types (including other vectors)


## R Variables are Objects

* The variables in R are technically known as objects
    - Objects should have meaningful names
    - Try to avoid common function names such as mean and sqrt or else it gets confusing
* Object names CANNOT start with a number
* Object names CAN have “.” and numbers within them
* Avoid “_”

**Variables (objects) can be created by = or <-**
```
x <- 'abc'
x = 'abc'
x
```

*abc*

**Determine data type of a variable**
```
typeof(x)
```

*character*


**Determine the number of variables in a variable**
```
length(x)
```

*1*

**Determine the number of characters in a variable**
```
nchar(x)
```

*3*


### R Data Types


* Character
    - “a”, “abc”,”hello world”
* Numeric (double or integer)
    - 1, 10.3, -199
* Logical
    - TRUE, FALSE
* Complex
    - 1+4i

**These data types allow you to know what sort of functions can be performed** 

For example numbers can be added but characters can’t be
```
1 + 2
```
*[1] 3*

```
'a' + 'b'
```
*Error in "a" + "b" : non-numeric argument to binary operator*

### Determine data types

```
x <- 4
typeof(x)
```
*double*

```
y <- 4.8 
typeof(y)
x+y
```

*double*

*8.8*

***These data types allow you to know what sort of functions can be performed***

```
x <- 'a'
typeof(x)
y <- 'Hello There'
typeof(y)
```
*character*

*character*

**For example numbers can be added but characters can’t be**

```
x+y
```

*Error in x + y: non-numeric argument to binary operator*

```
x <- TRUE
typeof(x)
```
*logical*

```
y <- FALSE
typeof(y)
```

*logical*

```
x <- 1+4i
typeof(x)
```

*complex*

```
x <- c(1:4)
typeof(x)
```
*integer*


### Basic Arithmetic

**Numbers in R need not be given object names**

```
20+3
20-3
20*3
```
*23*

*17*

*60*

```
20*3
20/3
20^3
```

*60*

*6.66666666666667*

*8000*

**Remainder of the Division**

```
20 %% 3
```

*2*

**Integer of the Division**

```
20 %/% 3
```

*6*

```
z <- 3*3
sqrt(z)
```

*3*

```
x <- 20
y <- 3
x + y
x - y
x * y
x / y
x ^ y
x %% y 
x %/% y #(integer of the division)
```
*23*

*17*

*60*

*6.66666666666667*

*8000*

*2*

*6*

## Built-In Functions

* Built-in functions are operations that one can “perform” on object that are available in R
* User-defined functions are functions that are written by the user
* Packages are R functions that are written by the R community that need to be loaded before using them

### Getting Help with Functions

* R Help
    - help()
    - ?
* The help() function and ? help operator in R provide access to the documentation pages for R functions, data sets, and other objects, both for packages in the standard R distribution and for contributed packages. 
* To access documentation for the standard lm (linear model) 
    - help(lm)
    - help(“lm")
    - ?lm
    - ?"lm" (i.e., the quotes are optional).
* To access help for a function in a package that’s not currently loaded, specify in addition the name of the package: for the rlm() (robust linear model) function in the MASS package:
    - help(rlm, package="MASS")

### Built-In Math Functions

| Function | Description |
| ----------- | ----------- |
| abs(x) | absolute value |
| sqrt(x) | square root |
| ceiling(x) | ceiling(3.475) is 4 |
| floor(x) | floor(3.475) is 3 |
| trunc(x) | trunc(5.99) is 5 |
| round(x, digits=n) | round(3.475, digits=2) is 3.48 |
| signif(x, digits=n) | signif(3.475, digits=2) is 3.5 |
| cos(x), sin(x), tan(x) | also acos(x), cosh(x), acosh(x), etc. |
| log(x) | natural logarithm |
| log10(x) | common logarithm |
| exp(x) | e^x |

**absolute value**

```
x <- -2
abs(x)
```
*2*

**square root**
```
x <- 4
sqrt(x)
```
*2*

**Rounding and Creating Integers**
```
x <- 4.693959
n <- 3
ceiling(x)
```
*5*
```
floor(x)
```
*4*
```
trunc(x)
```
*4*
```
round(x, digits=n)
```
*4.694*
```
signif(x, digits=n)
```
*4.69*

**Triganometry**
```
x <- 60
cos(x)
```
*-0.952412980415156*
```
sin(x)
```
*-0.304810621102217*
```
tan(x)
```
*0.320040389379563*

**Logarithms and Exponents**
```
x <- 3.2
log(x)
```
*1.16315080980568*
```
log10(x)
```
*0.505149978319906*
```
exp(x)
```
*24.5325301971094*


## Built In Statistical Functions

### Means, Medians, Ranges and Other Basic Functions on Number Sets


```
x <- 1:20
mean(x)
sd(x)
median(x)
min(x)
max(x)
```

```
quantile(x)
```

| 0% | 25% | 50% | 75% | 100% |
|------:|------:|------:|------:|------:|
| 1.00 | 5.75 | 10.50 | 15.25 | 20.00 | 

```
range(x)
sum(x)
diff(x, lag=1)
y <- scale(x, center=TRUE, scale=TRUE)
plot(x,y)
```
### Normal Distribution

* dnorm(x)
* pnorm(q)
* qnorm(p)
* rnorm(n, m=0,sd=1)

**normal density function (by default m=0 sd=1)**

```
x <- pretty(c(-3,3), 30)
y <- dnorm(x)
plot(x, y, type='l', xlab="Normal Deviate", ylab="Density", yaxs="i")
```

**cumulative normal probability for q**
```
pnorm(1.96)
```

**normal quantile: value at the p percentile of normal distribution **
```
qnorm(.9)
```

**n random normal deviates with mean m and standard deviation sd.  50 random normal variates with mean=50, sd=10 **
```
x <- rnorm(50, m=50, sd=10)
x
```

### Binomial Distribution

binomial distribution where size is the sample size and prob is the probability of a heads (for a coin toss) 

- dbinom(x, size, prob)
- pbinom(q, size, prob)
- qbinom(p, size, prob)
- rbinom(n, size, prob)

```
x <- seq(0,50,by=1)
y <- dbinom(x,50,0.2)
plot(x,y)
```

### Poisson Distribution

poisson distribution with m=std=lamda

- dpois(x, lamda)
- ppois(q, lamda)
- qpois(p, lamda)
- rpois(n, lamda)

```
x <- 0:20
y <- dpois( x=0:20, lambda=6 )
plot(x, y, xlim=c(-2,20))
```

### Uniform Distribution

- dunif(x, min=0, max=1)
- punif(q, min=0, max=1)
- qunif(p, min=0, max=1)
- runif(n, min=0, max=1)

```
numcases <- 10000
min <- 1
max <- 6
x <- as.integer(runif(numcases,min,max+1))
hist(x,main=paste( numcases," roles of a single die"),breaks=seq(min-.5,max+.5,1))

```

### Built-In String Functions

| Function | Description |
| ----------- |----------- |
| substr(x, start=n1, stop=n2) | Extract or replace substrings in a character vector. |
| grep(pattern, x , ignore.case=FALSE, fixed=FALSE) | Search for pattern in x. If fixed =FALSE then pattern is a regular expression. If fixed=TRUE then pattern is a text string. Returns matching indices. |
| sub(pattern, replacement, x, ignore.case =FALSE, fixed=FALSE) | Find pattern in x and replace with replacement text. If fixed=FALSE then pattern is a regular expression. |
| strsplit(x, split) | Split the elements of character vector x at split. |
| paste(..., sep="") | Concatenate strings after using sep string to seperate them |
| toupper(x) | Uppercase |
| tolower(x) | Lowercase |

### String Manipulation

**Concatenate strings after using sep string to seperate them.**
```
a <- "Hello"
b <- 'How'
c <- "are you? "
paste(a,b,c)
```
*[1] "Hello How are you? "*
```
paste(a,b,c, sep = "-")
```
*[1] "Hello-How-are you? "*
```
paste(a,b,c, sep = "", collapse = "")
```
*[1] "HelloHoware you? "*

```
paste('x',1:3,sep="")
```
*[1] "x1" "x2" "x3"*
```
paste('x',1:3,sep="M")
```
*[1] "xM1" "xM2" "xM3"*
```
paste('Today is', date())
```
*[1] "Today is Mon Dec 23 07:38:19 2019"*

**Text Capitalization: Create Upper and Lowercase strings**

```
a<- 'hEllo'
toupper(a)
```
*[1] "HELLO"*
```
tolower(a)
```
*[1] "hello"*

**Extract or replace substrings in a character vector.**

```
substring(a,1,2)
```
*[1] "He"*
```
substring(a,2,5)
```
*[1] "ello"*

```
x <- "abcdef" 
substr(x, 2, 4)
```
*bcd*

**grep**

*Search for pattern in x. If fixed =FALSE then pattern is a regular expression.  If fixed=TRUE then pattern is a text string. Returns matching indices.*

```
grep("A", c("b","A","c"), fixed=TRUE)
```
*2*

**substitution**

*Find pattern in x and replace with replacement text. If fixed=FALSE then pattern is a regular expression.  If fixed = T then pattern is a text string.*
```
sub("\\s",".","Hello There") 
```
*Hello.There*

**Split the elements of character vector x at split.**

```
strsplit("abc", "")
```
*[[1]]*<br>
*[1] "a" "b" "c"*

### Other Built in Functions

| Function | Description |
| ----------- | ----------- |
| seq(from , to, by) | generate a sequence |
| rep(x, ntimes | repeat x n times |
| cut(x, n) | cut divides the range of x into intervals and codes the values in x according to which interval they fall. The leftmost interval corresponds to level one, the next leftmost to level two and so on. | 	

**generate a sequence**
```
seq(1,10,2)
```
*[1] 1 3 5 7 9*

**Repeat a string**
```
rep(3,’a’)
```
*[1] "a" "a" "a"*

**divide continuous variable in factor with n levels**
```
# 
y <- cut(1:20, 20)
y
```

[1] (0.981,1.95] (1.95,2.9]   (2.9,3.85]   (3.85,4.8]   (4.8,5.75]   (5.75,6.7]   (6.7,7.65]   (7.65,8.6]  <br>
[9] (8.6,9.55]   (9.55,10.5]  (10.5,11.4]  (11.4,12.4]  (12.4,13.3]  (13.3,14.3]  (14.3,15.2]  (15.2,16.2] <br>
[17] (16.2,17.1]  (17.1,18.1]  (18.1,19.1]  (19.1,20]   <br>
20 Levels: (0.981,1.95] (1.95,2.9] (2.9,3.85] (3.85,4.8] (4.8,5.75] (5.75,6.7] (6.7,7.65] ... (19.1,20]


### Vectors and Matrices

A vector can be a “collection” of values or a single value

* Atomic vector
    - a collection of values
* Factors
    - special vectors that represent categorical data
* Matrix
    - special vector with rows and columns
* Data frame
    - a special data structure of rows and columns, the default structure for reading in “excel-like” files
* List
    - a vector of different data types (including other vectors)


**Vectors**

*An example of a numeric vector*

```
x <- 1
x
```
*[1]  1*

```
y <- c(1:10) 
y
length(y)
typeof(y)
```

*[1]  1  2  3  4  5  6  7  8  9 10*

*10*

*integer*

**Vectors can be any datatype (character, logical, complex)**

```
x
```
*[1] "a" "b" "c" "d"*

```
typeof(x)
```
*[1] "character"*

**For Numerical Vectors you can do comparisons**

```
x <- c(1:10)
x
```

*[1]  1  2  3  4  5  6  7  8  9 10*

```
x > 3
```
*[1] TRUE TRUE TRUE TRUE*

```
(x > 3) & (x < 8)
```
*[1] FALSE FALSE FALSE FALSE*

```
x[x > 3]
```
*[1] "a" "b" "c" "d"*

**Comparisons are logical objects**

```
typeof((x > 3) & (x < 8))
```

*logical*

**Combining Vectors**

```R
x <- c(1:10)
x <- c(x,11:20)
x
```
*[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20*

**Performing Math on Vectors**

```
x <- c(1:10)
x*3
```

* [1]  3  6  9 12 15 18 21 24 27 30*

```
x[4]
```
*[1] 4*

### Factors Are Nominal Vectors

The factor stores the nominal values as a vector of integers in the range [ 1... k ] and an internal vector of character strings (the original values) mapped to these integers.

```
genotype <- c(rep(“WT”,5),rep("KO",5))
factor(genotype)
```
*[1] WT WT WT WT WT KO KO KO KO KO*

*Levels: KO WT*

```
genotype <- factor(genotype,levels=c("WT","KO"))
genotype
```
*[1] WT WT WT WT WT KO KO KO KO KO*

*Levels: WT KO*

```
summary(genotype)
```
| WT | KO |
| ----------- | ----------- |
|5 | 5 |

## Matrix

**special vector with rows and columns**

* All columns in a matrix must have the same mode(numeric, character, etc.) and the same length. The general format is:
  - mymatrix <- matrix(vector, nrow=r, ncol=c, byrow=FALSE)
  - byrow=TRUE indicates that the matrix should be filled by rows. 
  - byrow=FALSE indicates that the matrix should be filled by columns (the default)


```
y <- matrix(1:20, nrow=5,ncol=4)
y 
cells <- c(1,26,24,68)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2") 
x <- matrix(cells, nrow=2, ncol=2, byrow=TRUE, dimnames=list(rnames, cnames))
x 
```

```R
x[,1] # 1st column of matrix
x[2,] # 2nd row of matrix 
x[2,1] # row 2, column 1
x[1,2] 
```

## Functions to Combine and Calculate Statistics on Matrices

* cbind(A,B,...)
    - Combine matrices(vectors) horizontally. Returns a matrix.
* rbind(A,B,...)
    - Combine matrices(vectors) vertically. Returns a matrix.
* rowMeans(A)
    - Returns vector of row means.
* rowSums(A)
    - Returns vector of row sums.
* colMeans(A)
    - Returns vector of column means.
* colSums(A)
    - Returns vector of column sums.
* t(A)
    - Transpose


```
y <- matrix(1:20, nrow=5,ncol=4,byrow = FALSE)
y
y <- matrix(1:20, nrow=5,ncol=4,byrow = TRUE)
y
t(y) #transpose
rowSums(y)
colMeans(y)
```

**Math with Matrices**

```
y*4
y*y
y/y
```

| Function | Description |
|------:|:-----|
| A * B	| Element-wise multiplication |
| A %*% B | Matrix multiplication |
| A %o% B | Outer product. AB' |
| crossprod(A,B) | A'B |
| crossprod(A) | A'A |
| t(A) | Transpose |
| diag(x) | Creates diagonal matrix with elements of x in the principal diagonal |
| diag(A) | Returns a vector containing the elements of the principal diagonal |
| diag(k) | If k is a scalar, this creates a k x k identity matrix. Go figure. |
| solve(A, b) |	Returns vector x in the equation b = Ax (i.e., A-1b) |
| solve(A) | Inverse of A where A is a square matrix. |

# Data Frames

* A special data structure of rows and columns, the default structure for reading in “excel-like” files
* A data frame is more general than a matrix, in that different columns can have different modes (numeric, character, factor, etc.). This is similar to SAS and SPSS datasets.

```
d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
x <- data.frame(d,e,f)
names(x) <- c("ID","Color","Passed")
x
```

```
setwd("~/Desktop/")
tbl <- read.csv(file="sample_example_R1_data_structures.csv",header=TRUE)
tbl
```

```
tbl[3:5] 
```

```
tbl[c("SampleID","Tissue")]
```

```
tbl$Gender 
```

```
tbl[tbl$SampleGroup == 'monocytes',]
```

```
subset(x=tbl,SampleGroup == 'monocytes',select=c('Tissue','SampleID'))
```

```
tbl1 <- read.csv(file="sample_example_R1_data_structures.csv",header=TRUE)
tbl2 <- read.csv(file="table2.csv",header=TRUE)
tbl2
```

## DataFrame Functions

* gives a very brief description of the data
    - str(df)
* gives the name of each variables
    - names(df)
* gives some very basic summary statistics for each variable
    - summary(df)
* shows the first few rows
    - head(df)
* shows the last few rows.
    - tail(df)	
* looks at duplicated elements and returns a logical vector. You can use table() to summarize this vector.
    - duplicated()
* keeps only the unique lines in a dataset
    - unique()	


```
merge.tbl <- merge(tbl1,tbl2,by='SampleID')
merge.tbl
```

```
help('merge')
```

## Lists

An ordered collection of objects (components). A list allows you to gather a variety of (possibly unrelated) objects under one name.


```
w <- list(name="Fred", mynumbers=a, mymatrix=y, age=5.3)
```

# Workspace Functions

* lists the objects in your workspace
    - ls()	
* removes an object in your workspace
    - rm(object1,object2)
* removes all objects in your workspace
    - rm(list=ls())
* saves R objects to a file
    - save(object1,object2,file=“file.RData”)
* load an R object from a file
    - load(“file.Rdata”)
* find current working directory   
    - getwd()
* set working directory    
    - setwd(‘C:/workingDirectory’) #	
* quit
    - quit()
* list all packages available to load
    - library()
* load package
    - library(package)
    - require(package)
* Install Packages
    - install.packages(“ggplot2")
    - source(“http://bioconductor.org/biocLite.R")
    - biocLite("DESeq2")
