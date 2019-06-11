# Main ACTIONet repository



## Installation of R package

### Linux

Install igraph c++ interface

```{bash}
sudo apt-get install libigraph0-dev
```

Or install it from the source from https://igraph.org/c/.

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```



### Mac OS

Upgrade LLVM to ensure OpenMP support:

```{bash}
brew install llvm
```

Install igraph c++ interface

```{bash}
brew install igraph
```

Or install it from the source from https://igraph.org/c/.

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```

