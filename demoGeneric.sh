 #!/bin/sh

TSNE_PERPLEXITY=6

if [ ! -f ./Generic.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/Generic.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf Generic.tar.gz
fi

mkdir -p output_generic
Rscript --vanilla co-fuse.R 'Generic' './Generic/GROUP1' './output_generic/group1' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R 'Generic' './Generic/GROUP2' './output_generic/group2' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'Generic' './Generic/GROUP1' './Generic/GROUP2' './output_generic/fisher_test'

