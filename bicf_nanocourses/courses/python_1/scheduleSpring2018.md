# [BICF](http://www.utsouthwestern.edu/labs/bioinformatics/) Python 1 - Nanocourse

This nanocourse will introduce Python for scientific computing.  
Python is an open-source, fun, easy to learn, and powerful programming language.
With deep community support and wide ranging deployment across many domains, Python is a worthy tool for projects large and small that any computational scientist should keep on hand.

Topics for this two day (all-day) course will include:

- Basic install, setup, and IDEs
- Basic Syntax
- Conditional statements, loops, functions
- Modules, classes, scripting, debugging
- Numerical arrays/matrices (numpy/scipy)
- Data structures (pandas)
- Plotting data (matplotlib, seaborn, bokeh)

The course will be interactive, with lectures followed by hands-on learning and exercises.
No previous programming experience is necessary. Familiarity with basic programming/scripting concepts is helpful. Students will also have the opportunity to share their own technical challenges and explore as a class how python can help.

### Contacts
* Course Coordinator: [Andrew Jamieson](mailto:andrew.jamieson@utsouthwestern.edu)
* Course Administration: [Rebekah Craig](mailto:rebekah.craig@utsouthwestern.edu)

### Preparation for Class

- Please bring your own laptop
- [Download & Install: Anaconda Python (with Python 3.6)](https://www.anaconda.com/download/)

	- Once installed, we will cover initial startup /setup of Python IDEs/ etc.. in the course.

- [Microsoft Teams - Python I Nanocourse Page](https://teams.microsoft.com/l/team/19%3a7010286a46384e80b00e912136e0f03d%40thread.skype/conversations?groupId=1056ac69-541c-4966-af52-1d4220a5d88e&tenantId=9d418695-71ac-4c31-b5b2-c196c8ec3c8a)

	- Log on to Teams to ask questions, share ideas, post content, get assignments
	- Use a browser (like Chrome) or download the desktop app


### Environment

Install [conda](https://conda.io) environment specified in [`environment.yml`](environment.yml) by downloading the file and running the following in Terminal (macOS and Linux) or an Anaconda Prompt (Windows):

```sh
# Install the environment
conda env create --file environment.yml
```

Activate the new environment (assumes `conda` version of [at least](https://github.com/conda/conda/blob/9d759d8edeb86569c25f6eb82053f09581013a2a/CHANGELOG.md#440-2017-12-20) 4.4):

* Windows: `activate bicf-python1`
* macOS and Linux: `source activate bicf-python1`

The environment should successfully install on both Linux, macOS and Windows.

Test the new environment by downloading the python script [`check_versions.py`](check_versions.py) and running:


```sh
# Test the environment
python check_versions.py
```
```
Output of script if everything works should be:

[  INFO] Hello - we're checking if your system is ready for the Python 1 Nanocourse
[  INFO] Python verion OK!
[  INFO] Checking for numpy 1.14.0
[  INFO] Checking for scipy 1.0.0
[  INFO] Checking for pandas 0.22.0
[  INFO] Checking for matplotlib 2.1.2
[  INFO] Checking for seaborn 0.8.1
[  INFO] Checking for bokeh 0.12.13
[  INFO] Checking for spyder 3.2.6

Woo! - Ready to go, see you at the nanocourse :-)
```

Deactivate the environment using:

* Windows: `deactivate bicf-python1`
* macOS and Linux: `source deactivate bicf-python1`

# Schedule

Day 1  | **February 27th, 2018**  
Room NB2.100A

| Time  | Topic | Instructor |
| ------------- | ------------- | ------------- |
| 9:00 - 10:00 a.m. | [Intro, IDEs, Setup](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=0.%2520Intro%2520%252B%2520Setup)  | Andrew Jamieson |
| 10:00 - 11:00 a.m.  | [Basic Syntax](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=1.%2520Basic%2520Syntax)  | Daniel Moser |
| 11:00 a.m. - 12:00 p.m.| [Practical Exercises](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=1.%2520Basic%2520Syntax) |  |
| 1:00 - 2:00 p.m.|[Control Statments, Loops, Functions](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=2.%2520Control-Functions-loops)| Daniel Moser |
| 2:00 - 3:00 p.m.|[Practical Exercises](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=2.%2520Control-Functions-loops) |  |
| 3:00 - 3:30 p.m.|[Modules, Classes, Environments](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=3.%2520Modules%2520%2526%2520Classes) | Benjamin Wakeland |
| 3:30 - 4:00 p.m.|[Scripting, Debugging](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=3.%2520Modules%2520%2526%2520Classes) | |
| 4:00 - 5:00 p.m.|[Practical Exercises](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=3.%2520Modules%2520%2526%2520Classes) | |

Day 2  | **February 28th, 2018**  
Room NB2.100A

| Time  | Topic | Instructor |
| ------------- | ------------- | ------------- |
| 9:00 - 10:00 a.m. | [Numpy + Scipy](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=4.%2520Numpy%2520%252B%2520SciPy)  | Viren Amin |
| 10:00 - 11:00 a.m.  | [Practical Exercises](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=4.%2520Numpy%2520%252B%2520SciPy)  |  |
| 11:00 a.m. - 12:00 p.m.| [Data Structures + Pandas](lectures/introduction_to_pandas_and_dataframes.slides.html) <br> [Slides and Workshop](lectures/pandas_workshop.zip)| Venkat Malladi |
| 1:00 - 2:00 p.m.|[Data Structures + Pandas (cont.)](lectures/introduction_to_pandas_and_dataframes.slides.html) <br> [Slides and Workshop](lectures/pandas_workshop.zip)| Venkat Malladi |
| 2:00 - 2:30 p.m.|[Plotting Data: Matplotlib, Notebook](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=6.%2520Matplotlib%2520%252B%2520Bokeh%2520%252B%2520Seaborn) | Wei Guo |
| 2:30 - 3:00 p.m.|[Plotting Data: Bokeh, Seaborn](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=6.%2520Matplotlib%2520%252B%2520Bokeh%2520%252B%2520Seaborn) |  |
| 3:00 - 4:00 p.m.|[Practical Exercises](https://teams.microsoft.com/_#/files/General?threadId=19%3A7010286a46384e80b00e912136e0f03d%40thread.skype&ctx=channel&context=6.%2520Matplotlib%2520%252B%2520Bokeh%2520%252B%2520Seaborn) |  |
| 4:00 - 5:00 p.m.|[Python Therapy: Student Cases](img/py.png) | All Instructors |

__TAs:__ Guillaume Jimenez, Gervaise Henry, Behrouz Saghafi Khadem

#### Python Therapy for Students

Opportunity to apply what you've learned to your own research! Students are encouraged to present a specific technical challenge encountered in their work and how they might solve this problem with python. We will review these cases as a class and discuss.

- [Student Data Upload Link - BioHPC Cloud](https://cloud.biohpc.swmed.edu/index.php/s/bdScPpoVlOsa6qx)
- [Teams Nanocourse - Python Therapy Channel - File Tab Upload](https://teams.microsoft.com/_#/files/Python%20Therapy?threadId=19%3Ade66d0e1eeb34888b1cd5797a8182d97%40thread.skype&ctx=channel)

## Additional Resources

- [PyCharm Community](https://www.jetbrains.com/pycharm/download/)
- [Offical Python Tutorial](https://docs.python.org/3/tutorial/)
- [Conda Cheatsheet](https://conda.io/docs/_downloads/conda-cheatsheet.pdf)
