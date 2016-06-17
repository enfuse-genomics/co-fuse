 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./DeFuse.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/DeFuse.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf DeFuse.tar.gz
fi

mkdir -p output_defuse
Rscript --vanilla co-fuse.R 'DeFuse' './DeFuse/GROUP1' './output_defuse/group1' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R 'DeFuse' './DeFuse/GROUP2' './output_defuse/group2' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'DeFuse' './DeFuse/GROUP1' './DeFuse/GROUP2' './output_defuse/fisher_test'

