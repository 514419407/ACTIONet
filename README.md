# Main ACTIONet repository



## Installation

### Linux

Install (upgrade to) the latest version of gcc compiler:

```{bash}
sudo add-apt-repository ppa:jonathonf/gcc-9.0
sudo apt-get install g++-9
```



Install igraph c++ interface

```{bash}
sudo apt-get install libigraph0-dev
```

Or install it from the source from https://igraph.org/c/.



Create  *~/.R* folder and copy *Makevars_gcc.txt* there under the name *Makevars*: 

```{bash}
 mkdir -p ~/.R
 wget --no-check-certificate https://raw.githubusercontent.com/shmohammadi86/ACTIONet/master/Makevars_gcc -O ~/.R/Makevars
```

Install Bioconductor dependencies:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install('scran', 'scater')
```

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```



### Mac OS

Install (upgrade to) the latest version of gcc compiler:

```{bash}
brew install gcc@9
```

Install igraph c++ interface

```{bash}
brew install igraph
```

Or install it from the source from https://igraph.org/c/.

Create  *~/.R* folder and copy *Makevars_gcc.txt* there under the name *Makevars*: 

```{bash}
 mkdir -p ~/.R
 wget --no-check-certificate https://raw.githubusercontent.com/shmohammadi86/ACTIONet/master/Makevars_gcc -O ~/.R/Makevars
```

Install Bioconductor dependencies:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install('scran', 'scater')
```

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```

