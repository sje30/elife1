# elife1

example one:  <https://elifesciences.org/articles/45952/figures>

Steps:


# Figure 1

data:

<https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDU5NTIvZWxpZmUtNDU5NTItZmlnMS1kYXRhMS12Mi5jc3Y=/elife-45952-fig1-data1-v2.csv?_hash=1t64b3sYektskcF74Wysrh9AgQcE6WMQq%2FsxKUrHnJk%3D>
code: <https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDU5NTIvZWxpZmUtNDU5NTItZmlnMS1jb2RlMS12Mi5y/elife-45952-fig1-code1-v2.r?_hash=QEAULzjYnZIQcxEFm%2Feoldo1PyS%2BX6iyz4GYXsTxCQM%3D>


Did not run:

```
source('elife-45952-fig1-code1-v2.r')
Error in eval(ei, envir) (from elife-45952-fig1-code1-v2.r#4) : object 'gtsB' not found
```


## Figure 2

Downloaded data and R

Again, didn't work.  Using an external package, `readxl` to read in
the data, but the file is not provided.  (See the file name: it is an
excel file, but the elife web page has an .csv file).  See also the
problem of the continuation character + appearing in the script.

```
source('elife-45952-fig2-code1-v2.r')
Error in source("elife-45952-fig2-code1-v2.r") : 
  elife-45952-fig2-code1-v2.r:4:24: unexpected '='
3: fig2c = read_excel("Fig 2C - YFP-fitness data.xlsx", 
4:             +     skip =
```


```
fig2c = read_excel("Fig 2C - YFP-fitness data.xlsx", 
            +     skip = 1)
```


## Figure 3

More woes.
```
source('elife-45952-fig3-code1-v2.r')
Error in read_excel("promoter-yfp data.xlsx", col_types = c("text", "numeric",  (from elife-45952-fig3-code1-v2.r#3) : 
  could not find function "read_excel"
```

Again, same problem, no .xlsx file provided and no readxl package
provided.

## Figure 4

Data files missing.  At least this one had some clear description of
packages that are required:

```
library(seqinr)
library(phangorn)
library(phytools)
```

But can see errors awaiting us... as none of the data flies were provided.


```
 9: my.tree = read.tree(file = "pruned tree")
12:my.alignment = read.alignment(file="gtsB alignment.fasta", format = "fasta")
77:mut.data = read.csv(file = "mutant fitness for comparison.csv")
