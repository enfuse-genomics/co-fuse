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

We provide the demo script which will download the data set and run our R scripts as described in the paper. Please run

```shell
./demo.sh
```

The downloaded data set consists of FusionCatcher results [1] on two different groups of samples: `AML` and `MM`. Each group consists of multiple samples. The first experiment (Recurrent Fusions) calls the R script `co-fuse.R` to compute recurrent fusions in both `AML` and `MM`. Experimental results can be found in the folder `./output_AML` and `./output_MM`.
The second experiment (Fisher's Exact test) calls the R scirpt `co-fuse2.R` to examine the relationship between gene fusions in `AML` and `MM`. Experimental results can be found in the folder `./output_fisher_test`.


## References

[1] Edgren,H., Murumagi,A., Kangaspeska,S., Nicorici,D., Hongisto,V., Kleivi,K., Rye,I.H., Nyberg,S., Wolf,M., Borresen-Dale,A.L. et al. (2011) Identification of fusion genes in breast cancer by paired-end RNA-sequencing. Genome Biol., 12, R6.

