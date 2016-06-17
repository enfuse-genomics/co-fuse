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

mkdir -p output_fusioncatcher
Rscript --vanilla co-fuse.R 'FusionCatcher' './FusionCatcher/AML' './output_fusioncatcher/AML' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R 'FusionCatcher' './FusionCatcher/MM' './output_fusioncatcher/MM' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'FusionCatcher' './FusionCatcher/AML' './FusionCatcher/MM' './output_fusioncatcher/fisher_test'
