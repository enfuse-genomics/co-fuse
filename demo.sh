 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./FusionCatcher.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/FusionCatcher.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf FusionCatcher.tar.gz
fi

Rscript --vanilla co-fuse.R './FusionCatcher/AML' './output_AML' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R './FusionCatcher/MM' './output_MM' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R './FusionCatcher' './output_fisher_test'

