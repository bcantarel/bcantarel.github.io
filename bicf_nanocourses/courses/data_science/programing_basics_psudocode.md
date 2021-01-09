| Concept                          | Example Pseudocode                         |
|----------------------------------|--------------------------------------------|
| Variables                        | x=3                                        |
|                                  | name="Bob"                                 |
|                                  | global userid = 123                        |
| Casting                          | str(3) returns "3"                         |
|                                  | int("3") returns 3                         |
|                                  | float("3.14") returns 3.14                 |
| Outputting to screen             | PRINT("hello")                             |
| Taking input from user           | name = INPUT("Please enter your name")     |
| Iteration – Count controlled     | FOR I = 0 to 7                             |
|                                  |    PRINT("Hello")                          |
|                                  | NEXT i                                     |
| Iteration – Condition controlled | WHILE answer != "computer”                 |
|                                  |    answer = INPUT("What is the password?") |
|                                  | ENDWHILE                                   |
|                                  | DO                                         |
|                                  |    Answer = INPUT("What is the password?") |
|                                  | UNTIL answer == "computer"                 |
| Logical operators                | WHILE x <=5 AND flag == FALSE              |
| Comparison operators             | ==                                         |
|                                  | !=                                         |
|                                  | <                                          |
|                                  | <=                                         |
|                                  | >                                          |
|                                  | >=                                         |
| Arithmetic operators             | +                                          |
|                                  | -                                          |
|                                  | *                                          |
|                                  | /                                          |
|                                  | MOD                                        |
|                                  | DIV                                        |
|                                  | ^                                          |
| Selection       | IF entry == "a" THEN  <br>
|                 |       PRINT("You selected A")  <br>
|                 | ELSEIF entry == "b" then   <br>                                |
|                 |    PRINT("You selected B")  <br>                               |
|                 | ELSE                                                       |
|                 |    PRINT("Unrecognised selection")                         |
|                 | ENDIF                                                      |
|                 | SWITCH ENTRY:                                              |
|                 |    CASE "A":                                               |
|                 |       PRINT("You selected A")                              |
|                 |    CASE "B":1                                              |
|                 |       PRINT("You selected B")                              |
|                 |    DEFAULT:                                                |
|                 |       PRINT("Unrecognised selection")                      |
|                 | ENDSWITCH                                                  |
| String handling | stringname.LENGTH                                          |
|                 | stringname.SUBSTRING(startingPosition, numberOfCharacters) |
| Subroutines     | FUNCTION triple(number)                                    |
|                 |    RETURN number * 3                                       |
|                 | ENDFUNCTION                                                |
|                 |                                                            |
|                 | Called from main program                                   |
|                 |                                                            |
|                 | Y =triple(7)                                               |
|                 | PROCEDURE greeting(name)                                   |
|                 |    PRINT("hello" + name)                                   |
|                 | ENDPROCEDURE                                               |
|                 |                                                            |
|                 | Called from main program                                   |
|                 |                                                            |
|                 | greeting("Hamish")                                         |
