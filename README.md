# pmsignature

The R package **pmsignature** is developed 
for efficiently extracting characteristic mutation patterns (mutation signatures) 
from the set of mutations collected typically from cancer genome sequencing data.

For extracting mutation signatures, 
principal component analysis or nonnegative matrix factorization have been popular.
Compared to these existing approaches, the **pmsignature** has following advantages:
  
  
1. **pmsignature** can perform robust estimation of mutation signatures in case of many contextual factors are taken into account such as two bases 5' and 3' to the mutated sites.
2. **pmsignature** provides intuitively interetable visualization of mutation signatures, which is analogous to sequencing logos.


## Paper

> Shiraishi et al. A simple model-based approach to inferring and visualizing cancer mutation signatures, PLoS Genetics, 2015,
[http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005657](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005657).

## Input data

Here, *mutation features* are elements used for categorizing mutations such as: 
  
* 6 substitutions (C>A, C>G, C>T, T>A, T>C and T>G)
* 5’ and 3’ flanking bases (A, C, G and T)
* transcription direction.

Currently, **pmsignature** can accept following two formats of tab-delimited text file.


### Mutation Position Format

sample1 chr1  100	A	C	
sample1	chr1	200	A	T	
sample1	chr2	100	G	T	
sample2	chr1	300	T	C	
sample3	chr3	400	T	C	
  
* The 1st column shows the name of samples 
* The 2nd column shows the name of chromosome 
* The 3rd column shows the coordinate in the chromosome
* The 4th column shows the reference base (A, C, G, or T).
* The 5th colum shows the alternate base (A, C, G, or T).


### Mutation Feature Vector Format

1 4	4	4	3	3	2	 
2	4	3	3	1	1	2	
3	4	4	3	2	2	2	
4	3	3	2	3	3	1	
5	3	4	2	4	4	2	
6	4	1	4	2	1	2	
3	2	1	1	1	1	2	
7	4	2	2	4	3	2	
  
* The 1st column shows the name of samples 
* From the 2nd to the last column show the value of mutation features.
* In this example, substitution patterns (1 to 6 values, C>A, C>G, C>T, T>A, T>C and T>G), two 5' and 3' flanking bases (1 to 4 values for each 4 site, A, C, G and T)
and transcription direction (1 to 2 values, + and -) are considered.
* However, you can put your favorite mutation features in this file format.

## Workflow

### Install the package

First, several R packages such as **ggplot2**, **Rcpp**, **GenomicRanges**, **BSgenome.Hsapiens.UCSC.hg19**,
which **pmsignature** depends has to be installed.
Also, **devtools** may be necessary for ease of installation.
Since the **pmsignature** utilizes C++ codes by way of **Rcpp**, you need to install C++ compiler
(e.g., [Rtools](http://cran.r-project.org/bin/windows/Rtools/) for Windows, Xcode for Mac).
See [Advanced R](http://adv-r.had.co.nz/Rcpp.html) by Dr. Hadley Wickham about the basic usage of **Rcpp**.
```
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19"))
install.packages("devtools")
install.packages("ggplot2")
install.packages("Rcpp")
```

Currently, the easiest way for installing **pmsignature** is to use the package **devtools**:
  
```
library(devtools)
devtools::install_github("friend1ws/pmsignature")
library(pmsignature)
```

For those who failed to installing, we recommend to use the newest R version.
Also, you may be required to upgrade the bioconductor:

```
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
```


### Prepare input data

First, create the input data from your mutation data.

After installing **pmsignature**,
you can find example files at the directory where pmsignature is installed:

* Mutation Position Format
```
inputFile <- system.file("extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz", package="pmsignature")
print(inputFile)
```

* Mutation Feature Vector Format
```
inputFile <- system.file("extdata/Hoang_MFVF.ind.txt.gz", package="pmsignature")
print(inputFile)
```


### Read input data

Type the following commands:
  
* Mutation Position Format
```
G <- readMPFile(inputFile, numBases = 5)
```
Here, *inputFile* is the path for the input file.
*numBases* is the number of flanking bases to consider including the central base (if you want to consider two 5' and 3' bases, then set 5).
You can format the data as the full model by typing
```
G <- readMPFile(inputFile, numBases = 5, type = "full")
```
Also, you can add transcription direction information by typing (in that case, the package **TxDb.Hsapiens.UCSC.hg19.knownGene** is necessary)
```
G <- readMPFile(inputFile, numBases = 5, trDir = TRUE)
```

* Mutation Feature Vector Format
```
G <- readMFVFile(inputFile, numBases = 5, type="independent", trDir=TRUE)
```

### Estimate the parameters


When you want to set the number of mutation signature as 3, type the following command (see also ?getPMSignature):
  
```
Param <- getPMSignature(G, K = 3)
```

If you want to add the background signature, then after obtaining the background probability, perform the estimation.
Currently, we only provide the background data for the "independent" and "full" model with 3, 5, 7 and 9 flanking bases (see also ?readBGFile).

```
BG_prob <- readBGFile(G)
Param <- getPMSignature(G, K = 3, BG = BG_prob)
```

In default, we repeat the estimation 10 times by changing the initial value,
and select the parameter with the highest value of log-likelihood.
If you want to changet the trial number, then

```
Param <- getPMSignature(G, K = 3, numInit=20)
```


### Visualing the mutation signatures and memberships

To visualize the mutation signatures by typing (see also ?visPMSignature):
```
visPMSignature(Param, 1)
visPMSignature(Param, 2)
visPMSignature(Param, 3)
```

To obtain the value of estimated mutation signatures (see also ?getSignatureValue):
```
getSignatureValue(Param, 1)
getSignatureValue(Param, 2)
getSignatureValue(Param, 3)
```


To see the overview of the estimated membership parameter (see also ?visMembership):
```
visMembership(G, Param)
```
Here, samples are sorted according to the number of mutations.
To unsort the sample, set sortSampleNum = FALSE.
Also, not to multiply the number of mutations to the barplot,
set multiplySampleNum = FALSE. 

To obtain the value of estimated membership parameters (see also ?getMembershipValue):
```
getMembershipValue(Param)
```
