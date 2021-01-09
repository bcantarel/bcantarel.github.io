## R Nanocourse - Programming basics

### Functions, conditional expressions, iterations


#### Functions

If you need to perform the same operations over and over again.
You will encounter situations in which the function does not already exist, so language permit you to write your own.

Basics of a function:

1. A function is a named sequence of statements that performs a computation.
2. A function may or may not accept an argument(s).
3. A function may or may not return an object.

Function Syntax

```
Function name
    Pass In: variables
    Do something
    Pass Out: value
Endfunction
```

Let's start with some easy examples:

**Question 1: Write a function to find the square of any number?**

**Question 2: Write a function converts temperatures from Fahrenheit to Celsius?**



### Conditional expressions

Conditional expressions are one of the basic features of programming.
They are used for what is called *flow control*.
The most common conditional expression is the if-else statement.

#### Boolean expressions

1. A boolean expression is an expression that is either true or false.
2. Boolean expressions allow the programmer to create controls structures.
3. Expressions can be tested to assess whether they are true or false

```
a = 5
b = 5
```

**Question 3: Is a equal to b?**

```
a == b
```

**Question 4: Is a greater than b?**

```
a > b
```

#### Logical operators

Logical operators are used with relational operators to create more complex boolean expressions.

* and
* or
* not

```
a = 100
b = 10 x 9

c = 4
d = 2 + 2
```

**Question 5: Is a equal to b?**

**Question 6: Is c equal to d?**

**Question 7: Is a equal to b and c equal to d?**

```
True if both expressions are True

a == b and c == d
```

**Question 8: Is a equal to b or c greater than d?**


```
True if one expression is True

a == b or c > d
```

**Question 9: Is a not equal to d?**

```
True if the boolean expression is False, or False if the boolean expression is True

not (c == d)
```

#### Conditional execution
* Boolean expressions by themselves are not that useful.
* Use the results of boolean expressions to control the execution of a program


#### IF Statement
```
IF something
    THEN do something
ENDIF
```

**Question 10: If a number is odd square the number?**


#### Alternative conditional execution
* A second form of the if statement is “alternative execution”, in which there are two possibilities, True or False, and the condition determines which one runs.

#### IF-ELSE Statement
```
IF something
    THEN do something
ENDIF
```

**Question 11: If a number is odd, square the number or divide by 2?**


#### Chained conditionals
* It may be useful to consider multiple conditionals.

```
IF something
  THEN do something
ELSE IF sommething else
    then do random thing
ELSE
    do another thing
ENDIF
```

**Question 11: You have test scores from 0-100. Give a letter grade based on where a students score falls?**


#### Nested conditionals
* A conditional can also be nested within another.

```
IF something
  THEN do something
ELSE
    IF something b
      THEN do this
ENDIF
```

**Question 12: You have test scores from 0-100. Give a letter grade based on where a students score falls.
If the student gets below 50, tell let the teacher know that this student will fail.**

#### Iterations
* Computers are often used to automate repetitive tasks.
* Repeating identical or similar tasks without making errors is something that computers do well and people do poorly.
* In a computer program, repetition is also called iteration.


#### For statment
* The For loop takes a group of elements and runs the code within the loop for each element.

```
FOR condition
  do something
ENDFOR
```

**Question 13: For all test scores in a class of 10, give a letter grade**

```
[80, 90, 30, 60, 99, 100, 56, 78, 67, 88]
```

**Question 14: For all test scores in a class of 10, give a letter grade. If the student gets below 50, tell let the teacher know that this student will fail.**


#### White statement
* The while statement executes a body of code while a condition is true

More formally, here is the flow of execution for a while statement:

1. Determine whether the condition is true or false.
2. If false, exit the while statement and continue execution at the next statement.
3. If the condition is true, run the body and then go back to step 1.

```
PRECONDITION: variable X is equal to 1
WHILE x < 5
  Compute x + 1
ENDWHILE
```

**Question 15: While a number is less than 100, square the number. **


#### Now lets combine everything

**Question 16: Make a function that computes the average grade for odd number test scores. **

```
[80, 90, 30, 60, 99, 100, 56, 78, 67, 88]
```
