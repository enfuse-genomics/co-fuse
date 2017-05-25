# Co-fuse: A recurrent fusions identification and analysis tool based on RNA-sequencing

Co-fuse implements a technique described in "Co-fuse: a recurrent fusions identification and analysis tool based on RNA-sequencing". 
Co-fuse identifies and compares recurrent fusion genes across different samples. It relies on the results of existing fusion gene detection tools. The program is written in R.

## Requirements

In order to use Co-fuse, RStudio needs to be installed.
RStudio can be installed by following this [link](https://www.rstudio.com/products/rstudio/download/).

Please follow this [link](http://www.r-bloggers.com/installing-r-packages/) for instuctions on how to install the following R Libraries/packages.


```
install.packages("dplyr")
install.packages("tsne")
install.packages("ggplots")
install.packages("rpart")
```

## Usage

### Recurrent Fusions

The algorithm identify a set of recurrent fusion genes occurred across multiple samples. Assuming that the output of fusion software can be found in the folder `SOFTWARE_RESULTS/muliple_samples` (e.g., TopHat/samples), the re-current fusion algorithm can be executed as follow:

```
Rscript --vanilla co-fuse.R 'software' 'SOFTWARE_RESULTS/muliple_samples' './output_dir' _TSNE_PERPLEXITY_
```

Here `'software'` is set to one of the following softwares: '_defuse_', '_fusioncatcher_', '_tophat_', '_soapfuse_' or '_generic_'.
`_TSNE_PERPLEXITY_` represents the perplexity parameter. It can be thought as a parameter that sets the number of effective nearest neighbors. A larger or denser dataset requires a larger perplexity. Typical values for the perplexity range between 5 and 50. We use the perplexity of 6 in our experiments.

The output will be stored in the folder `./output_dir`. It consists of multiple files:

- **summary.AB.csv** A CSV file listing recurrent genes and samples containing these recurrent genes
- **barplot.AB.pdf** Display recurrent genes sorted by their frequency 
- **heatmap.pdf** Cluster analysis to discover relationship between samples
- **tsne_plot.pdf** Dimensionality reduction for visualizing samples using t-Distributed Stochastic Neighbor Embedding (t-SNE)

### Fisher's Exact Test

Fisher’s exact test examines the relationship between the two dimensions of the contingency table. The null hypothesis is that the relative proportions of gene fusions in group 1 are the same as the relative proportions of gene fusions in group 2. In other words, these two groups are not different. 

For each fusion, we calculate a two-sided p-value. 
If the p-value (prbability) is small, the null hypothesis can be rejected. In other words, there is a significant
difference between two groups. The Fisher's exact test can be executed as follow:

```
Rscript --vanilla co-fuse2.R 'software' 'SOFTWARE_RESULTS/GROUP_1' 'SOFTWARE_RESULTS/GROUP_2' './output_dir'
```

Here `'software'` is set to one of the following softwares: '_defuse_', '_fusioncatcher_', '_tophat_', '_soapfuse_' or '_generic_'.
'SOFTWARE_RESULTS/GROUP_1' represents the directory containing fusion results from the first group.
'SOFTWARE_RESULTS/GROUP_2' represents the directory containing fusion results from the second group.
The output table (**twoGroups.csv** --- in the csv format) will be stored in the folder `./output_dir`.



## Demonstration

### Reproducing experimental results reported in the paper

We provide the demo script which will reproduce the results reported in our paper. The script will download the data set and run both R scripts as described in the paper. Please execute:

```shell
./demoFusionCatcher.sh
```

The downloaded data set contains FusionCatcher results (Nicorici et al., 2014) on two groups of samples: `AML` and `MM`. Each group consists of multiple samples. The first experiment (Recurrent Fusions) calls the R script `co-fuse.R` to compute recurrent fusions from both `AML` and `MM`. Experimental results can be found in the folder `./output_reproduce_results/AML_MM`.

The second experiment (Fisher's Exact test) tests whether there exists a significant difference between two groups of samples. It calls the R script `co-fuse2.R` to examine the relationship between gene fusions in `AML` and `MM`.
Experimental results (**twoGroups.csv**) will be stored in the folder `./output_reproduce_results/fisher_test`.


### FusionCatcher, SoapFuse, TopHat and DeFuse

We provide the demo script to illustrate the use of our software with four of popular fusion softwares: FusionCatcher (Nicorici et al., 2014), SoapFuse (Jia et al., 2013), TopHat (Kim and Salzberg, 2011) and DeFuse (McPherson et al., 2011). The demo script can be executed by:

```shell
./demoFusionCatcher.sh
./demoSoapFuse.sh
./demoTopHat.sh
./demoDeFuse.sh
```

### Generic fusion softwares

For other popular fusion softwares, we expect the first two columns to be Gene1 and Gene2. The demo script can be executed by:

```shell
./demoGeneric.sh
```

If you have any questions on use fusion softwares not described here, please feel free to contact us.


<!--
#### Note
1. Here we use FusionCatcher as our software. Other fusion algorithms can also be applied here. However, the code needs to be slightly modified. For example, line 59 in `co-fuse.R` would need to be modified: 

    ```
    filename <- list.files(path=folders[i],pattern="*.GRCh37.txt",full.name=TRUE,recursive=TRUE)
    ```

    Instead of setting `pattern="*.GRCh37.txt"`, one should set `pattern` to point to the output of other fusion softwares.

2. We assume that the first two columns contain a pair of recurrent genes. FusionCatcher output a pair of gene (geneA and geneB) in the first two columns. For other fusion software, one might need to modify the following R code:

    ```
    df <- data.frame(geneA=dat[,1],geneB=dat[,2],stringsAsFactors = F)
    ```
-->





### Known errors

```
./demoXXXXXX.sh: line XX: Rscript: command not found
```
Rstudio has not been installed. Please follow the instruction above to install RStudio.

## References

Jia, W., Qiu, K., He, M., Song, P., Zhou, Q., Zhou, F., et al. (2013). SOAPfuse:
an algorithm for identifying fusion transcripts from paired-end RNA-Seq data.
Genome biology, 14(2).

Kim, D. and Salzberg, S. L. (2011). TopHat-Fusion: an algorithm for discovery of
novel fusion transcripts. Genome biology, 12(8), R72.

McPherson, A., Hormozdiari, F., Zayed, A., Giuliany, R., Ha, G., Sun, M. G.,
Griffith, M., Moussavi, A. H., Senz, J., Melnyk, N., et al. (2011). deFuse: an
algorithm for gene fusion discovery in tumor RNA-Seq data. PLoS Computational
Biology, 7(5).


Nicorici, D., Satalan, M., Edgren, H., Kangaspeska, S., Murumagi, A., Kallioniemi,
O., et al. (2014). FusionCatcher - a tool for finding somatic fusion genes in paired-
end RNA-sequencing data. bioRxiv.

<!--
Agrawal, R., Imielinski, T., and Swami, A. (1993). Mining association rules between
sets of items in large databases. In Proc. ACM SIGMOD International Conference
on Management of Data, pages 207–216.

Barretina, J., Caponigro, G., Stransky, N., Venkatesan, K., Margolin, A. A., Kim,
S., et al. (2012). The cancer cell line encyclopedia enables predictive modelling
of anticancer drug sensitivity. Nature, 483(7391), 603–607.

Beccuti, M., Carrara, M., Cordero, F., Lazzarato, F., Donatelli, S., Nadalin, F.,
Policriti, A., and Calogero, R. A. (2014). Chimera: a bioconductor package for
secondary analysis of fusion products. Bioinformatics, 30(24), 3556–3557.

Capdeville, R., Buchdunger, E., .Zimmermann, J., and Matter, A. (2002). Glivec
(STI571, imatinib), a rationally developed, targeted anticancer drug. Nature
reviews Drug discovery, 1(7), 493–502.

der Maaten, L. V. and Hinton, G. (2008). Visualizing data using t-SNE. Journal of
Machine Learning Research, 9(85), 2579–2605.

Hoogstrate, Y., Bottcher, R., S, S. H., van der Spek, P. J., Jenster, G., and Stubbs, A. P.
(2016). FuMa: reporting overlap in rna-seq detected fusion genes. Bioinformatics,
32(8), 1226–1228.

Jia, W., Qiu, K., He, M., Song, P., Zhou, Q., Zhou, F., et al. (2013). SOAPfuse:
an algorithm for identifying fusion transcripts from paired-end RNA-Seq data.
Genome biology, 14(2).

Kim, D. and Salzberg, S. L. (2011). TopHat-Fusion: an algorithm for discovery of
novel fusion transcripts. Genome biology, 12(8), R72.

Latysheva, N. S. and Babu, M. M. (2016). Discovering and understanding oncogenic
gene fusions through data intensive computational approaches. Nucleic acids
research, 44(10), 4487–4503.

McPherson, A., Hormozdiari, F., Zayed, A., Giuliany, R., Ha, G., Sun, M. G.,
Griffith, M., Moussavi, A. H., Senz, J., Melnyk, N., et al. (2011). deFuse: an
algorithm for gene fusion discovery in tumor RNA-Seq data. PLoS Computational
Biology, 7(5).

Mertens, F., Johansson, B., Fioretos, T., and Mitelman, F. (2015). The emerging
complexity of gene fusions in cancer. Nature reviews Cancer, 15(6), 371–381.

Morgan, G. J., Walker, B. A., and Davies, F. E. (2012). The genetic architecture of
multiple myeloma. Nature reviews Cancer, 12(5), 335–348.

Nicorici, D., Satalan, M., Edgren, H., Kangaspeska, S., Murumagi, A., Kallioniemi,
O., et al. (2014). FusionCatcher - a tool for finding somatic fusion genes in paired-
end RNA-sequencing data. bioRxiv.

Shugay, M., de Mendibil, I. O., Vizmanos, J. L., and Novo, F. J. (2013). Oncofuse:
a computational framework for the prediction of the oncogenic potential of gene
fusions. Bioinformatics, 29(20), 2539–2546.

Wang, Q., Xia, J., Jia, P., Pao, W., and Zhao, Z. (2013). Application of next generation
sequencing to human gene fusion detection: computational tools, features and
perspectives. Briefings in bioinformatics, 14(4), 506–519.

Wilks, C., Cline, M. S., Weiler, E., Diehkans, M., Craft, B., Martin, C., et al. (2014).
The Cancer Genomics Hub (CGHub): overcoming cancer through the power of
torrential data. Database : the journal of biological databases and curation.

Wong, S. and Witte, O. N. (2004). The bcr-abl story: bench to bedside and back.
Annual review of immunology, 22, 247–306.
-->
