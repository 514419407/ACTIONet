# Main ACTIONet repository



### Installation

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```

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



Create  *~/.R* folder and copy *Makervars_gcc.txt* there under the name *Makevars*: 



```{bash}

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







### 