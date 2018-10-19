## R Nanocourse - Introduction to R for Beginners

#### Introduction to R Workshop

##### First we will need to set up your working directory

1.  Method 1
    1.  Setup a new folder and name it "RNanoCourse" on your desktop.
    2.  Download "Statistical Test Dataset", in your download directory, you will see dig_csv.zip, unzip it and copy DIG_csv folder into "RNanoCourse"
    3.  Start R studio, click File-New File-R Script and save the new R Script to the directory Desktop-RNanoCourse-DIG_csv.
    4.  Click Session-Set Working Directory-To Source File Location
2.  Method 2, you can save your R script anywhere you like
    1.  Setup a new folder and name it "RNanoCourse" on your desktop.
    2.  Download "Statistical Test Dataset", in your download directory, you will see dig_csv.zip, unzip it and copy DIG_csv folder into "RNanoCourse"
    3.  Click Session- Choose Directory - Desktop-RNanoCourse-DIG_csv
3.  Method 3 in R Studio
    *   setwd("~/Desktop/RNanoCourse/DIG_csv") #on Mac
    *   setwd("C:/Users/username/Desktop/RNanoCourse/DIG_csv") #on PC

##### Examples of functions that do not need any arguments

```
ls()  
library()  
getwd()  
colors()
```

##### Let's use R as a calculator

```
a<-2  
a  
b<-3  
a+b  
a^b  
a.vec<-c(2,4,6)  
a.vec  
length(a.vec)  
a.vec[1:2]
```

##### Creating a vector

```
b.vec<-c(3,5,7)  
a.vec+b.vec  
a.vec*b.vec  
A.vec
```

You will see an error as R is case sensitive

##### Creating a matrix

```
a.mat<-matrix(data=c(1,2,3,4,5,6),nrow=3,ncol=2)  
a.mat  
a.mat^2  
log2(a.mat)  
a.mat[1:2,1:2]  
a.mat[1:2,1]  
a.mat[1:2,]
```

##### Variables can contain strings

```
message<-"Hello world"  
message  
c.vec<-c('Hello','Goodbye')  
c.vec
```

##### The help function is useful to learn how to use other functions

To get to a function's help file, use ?function_name

```
?t.test  
?mean
```

Press q to return

##### Mathematical functions

```
x<-2.1  
exp(x)  
log(x,base=exp(1))  
log2(x)  
log10(x)  
sqrt(x)  
abs(x)  
round(x)  
factorial(x)  
choose(5,3)
```


```
cos(x)  
sin(x)  
tan(x)
```

##### Object types

```
x<-c(0,1,0,1)  
logical(length=3)  
is.logical(x)  
as.logical(x)  
is.character(x)  
as.character(x)  
is.numeric(x)  
as.numeric(x)  
is.integer(x)  
as.integer(x)
```

##### Working with vectors

```
my.vec <- c('h','e','l','l','o')  
my.vec  
my.vec1 <- seq(from=1,to=10,by=2)  
my.vec1  
my.vec2 <- 1:4  
my.vec2  
my.vec2 <- c(my.vec2, NA)  
my.vec2  
my.vec1[my.vec1<6]  
subset(my.vec1,my.vec1<6)  
plot(my.vec2,my.vec1)
```

##### Functions to learn about an object

```
length(my.vec)  
nchar('abcdefg')  
dim(a.mat)  
mode(my.vec)  
typeof(my.vec)  
attributes(my.vec)  
Res<-t.test(rnorm(20,1),rnorm(20,2))  
Res  
attributes(Res)  
Res$p.value  
length(my.vec)
```

##### Working with matrices

```
my.mat1<-matrix(1:20,nrow=5,ncol=4)  
my.mat1[!is.na(my.vec2),]  
my.mat2<- cbind(1:5,seq(1,10,2))  
my.mat2  
my.mat2<- cbind(a=1:5,b=seq(1,10,2))  
my.mat2
```

##### Working with lists

```
my.list <- list(fruits=c('apple','grape','banana'), num=3,colors=c('red','green','yellow'))  
my.list$fruits  
my.list[[1]]  
my.list[[1]][2]  
names(my.list)  
unlist(my.list)  
length(my.list)  
attributes(my.list)
```

##### Working with data frames

```
numbers <- 1:4  
letters<-c('a','b','c','d')  
grp<- as.factor(c(1,0,0,1))  
mydata <- data.frame(cbind(numbers, letters, grp))  
mydata # View the data frame  
mydata[1:2,]  
mode(mydata)  
attributes(mydata)  
colnames(mydata)
rownames(mydata)
```

##### Reading in a file

```
x <- read.table(file="dig.csv", sep=",", header=T)
head(x)
```

##### Installing a package

```
source("http://bioconductor.org/biocLite.R")  
biocLite("DESeq2")
install.packages("ggplot2")
```

##### Exiting R

```
q() # save workspace? (y/n)
```