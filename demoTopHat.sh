 #!/bin/sh

TSNE_PERPLEXITY=6

if [ ! -f ./TopHat.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/TopHat.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf TopHat.tar.gz
fi

mkdir -p output_tophat
Rscript --vanilla co-fuse.R 'TopHat' './TopHat/GROUP1' './output_tophat/group1' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R 'TopHat' './TopHat/GROUP2' './output_tophat/group2' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'TopHat' './TopHat/GROUP1' './TopHat/GROUP2' './output_tophat/fisher_test'

