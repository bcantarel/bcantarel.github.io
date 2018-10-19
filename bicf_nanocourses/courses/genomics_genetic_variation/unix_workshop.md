# Unix Basics Review

### Common Commands
| Command| Function |
| ------------- |:-------------:|
| ls | list contents|
| cd | change directory|
| mkdir | make a directory|
| rm | use caution, it is easy to delete more that you would like|
| head | prints the top few lines to the terminal window|
| tail | prints the last few lines to the terminal window|
| sort | sorts the lines|
| uniq | prints the unique lines|
| grep | filnds the lines that contain a pattern|
| wc | counts the number of lines, characters and words|
| mv | move files|
| cp | copy files|
| date | returns the current date and time|
| pwd | return working directory name|
| ssh | remote login|
| scp | remote secure copy|
| ~ | represents your home directory|
| man [command] | manual page for the command|

### text editors
1. emacs
2. vi
3. vim
4. nedit
5. nano

### Unix Problem Set
* Log into your machine or account.
  1. Mac Users: Open the application: Terminal
  2. From Terminal: `ssh -Y  student01@nucleus.biohpc.swmed.edu
  3. Windows users -- Please refer to the [PuTTY instructions]() with your username and the server `nucleus.biohpc.swmed.edu`
* What is the full path to your home directory?
* Go up one directory?
  1. How many files does it contain?
  2. How many directories?

* Using a text editor create a fasta file and name it sequences.fasta. Make sure it ends up in the proper directory, locally or remotely.
  1. This is fasta file format:

`>seqName description
ATGGCGTCTTGGCCTTAAAAGCTC`

* Without using a text editor examine the contents of the file sequences.fasta.
  1. How many lines does this file contain?
  2. How many characters?    (Hint: check out the options of wc)
  3. What is the first line of this file?    (Hint: read the man page of head)
  4. What are the last 3 lines?    (Hint: read the man page of tail)
  5. How many sequences are in the file?    (Hint: use grep)

* Rename sequences.fasta to something more informative of the sequences the file contains. (Hint: read the man page for mv)
* Create a directory called fasta.     (Hint: use mkdir)
* Copy the fasta file that you renamed to the fasta directory. (Hint: use cp)
* Verify that the file is within the fasta directory.    (Hint: use ls fasta/)
* Delete the the original file that you used for copying.    (Hint: use rm, be careful)
* Read the man page for rm and cp to find out how to remove and copy a directory.
* Print out your history and redirect it to a file called unixBasics.history.txt

### Commands to try

`ls -l
ls -lt`

You can string more than one command together with a pipe (|) , such that the output of the first command is received by the second command.

`ls -lt | head`

You can string more than one command together with a semi-colon (;) , such that the commands run sequentially, but that output does not get passed into the next command.

`date; some program command ; date`

You can redirect the output of a command into a file

`grep PATTERN > PATTERN.txt`

You can append the output of a command to a file

`grep PATTERN2 >> PATTERN.txt`

You can redirect stderr to a file

`command 2> filename`

You can redirect the output (stdout) and stderr to a file

`command &> filename`
