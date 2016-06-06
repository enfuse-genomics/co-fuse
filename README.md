# Co-fuse: A recurrent fusions identification and analysis tool based on RNA-sequencing

Co-fuse implements a technique described in Co-fuse: a recurrent fusions identification and analysis tool based on RNA-sequencing. The program is written in R.

## Requirements

In order to use Co-fuse, RStudio needs to be installed.
RStudio can be installed by following this [link](https://www.rstudio.com/products/rstudio/download/).

Please follow this [link](http://www.r-bloggers.com/installing-r-packages/) for instuctions on how to install the following R Libraries/packages.


```
install.packages("dplyr")
install.packages("tsne")
install.packages("ggplots")
```

## Demo

We provide the demo script which will download the data set and run our R scripts as described in the paper. Please execute

```shell
./demo.sh
```

The downloaded data set consists of FusionCatcher results [1] on two different groups of samples: `AML` and `MM`. Each group consists of multiple samples. The first experiment (Recurrent Fusions) calls the R script `co-fuse.R` to compute recurrent fusions in both `AML` and `MM`. Experimental results can be found in the folder `./output_AML` and `./output_MM`.
The second experiment (Fisher's Exact test) calls the R scirpt `co-fuse2.R` to examine the relationship between gene fusions in `AML` and `MM`. Experimental results can be found in the folder `./output_fisher_test`.

## Usage

### Recurrent Fusions

The algorithm identify a set of recurrent fusion genes occurred across multiple samples. Assuming that the output of FusionCatchers can be found in the folder `FusionCatcher/group`, the re-current fusion algorithm can be executed as follow.

```
Rscript --vanilla co-fuse.R './FusionCatcher/group' './output_group' _TSNE_PERPLEXITY_
```

Here `_TSNE_PERPLEXITY_` represents the perplexity parameter. It can be thought as a parameter that sets the number of effective nearest neighbors. A larger or denser dataset requires a larger perplexity. Typical values for the perplexity range between 5 and 50. We use the perplexity of 5 in our experiments.

The output will be stored in the folder `./output_group`. It consists of multiple files:
- **summary.AB.csv** A CSV file listing recurrent genes and samples containing these recurrent genes
- **barplot.AB.pdf** Display recurrent genes sorted by their frequency 
- **heatmap.pdf** Cluster analysis to discover relationship between samples
- **tsne_plot.pdf** Dimensionality reduction for visualizing samples using t-Distributed Stochastic Neighbor Embedding (t-SNE)

#### Note
1. Here we use FusionCatcher as our software. Other fusion algorithms can also be applied here. However, the code needs to be slightly modified. For example, line 59 in `co-fuse.R` would need to be modified:

```
filename <- list.files(path=folders[i],pattern="*.GRCh37.txt",full.name=TRUE,recursive=TRUE)
```

Instead of setting `pattern="*.GRCh37.txt"`, one should set `pattern` to point to the output of other fusion softwares.

2. We assume the first two columns contain a pair of recurrent genes. FusionCatcher output a pair of gene (geneA and geneB) in the first two columns. For other fusion software, one might need to modify the following R code:

```
df <- data.frame(geneA=dat[,1],geneB=dat[,2],stringsAsFactors = F)
```


### Fisher's Exact test

Fisher's exact test examines the relationship between the two dimensions of the contingency table. The algorithm tests if there exists a significant difference between two groups of samples.

```
Rscript --vanilla co-fuse2.R './FusionCatcher/group1' './FusionCatcher/group2' './output_fisher_test'
```

The output will be stored in the folder `./output_fisher_test`. It consists of:
- **twoGroups.csv**


<!---
TODO:
1. Fix co-fuse2.R to accept 3 arguments and check if the first 2 arguments contain more than one sample
2. Write more description on Fisher's Exact test
-->


## References

[1] Edgren,H., Murumagi,A., Kangaspeska,S., Nicorici,D., Hongisto,V., Kleivi,K., Rye,I.H., Nyberg,S., Wolf,M., Borresen-Dale,A.L. et al. (2011) Identification of fusion genes in breast cancer by paired-end RNA-sequencing. Genome Biol., 12, R6.

