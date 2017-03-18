 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./DeFuse.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://www.dropbox.com/s/d9g0xm6hgd3m1jf/DeFuse.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf DeFuse.tar.gz
    
    #
    # Merge GROUP 1 and GROUP 2
    #
    mkdir -p DeFuse/GROUP12
    cp -r DeFuse/GROUP1/* DeFuse/GROUP12
    cp -r DeFuse/GROUP2/* DeFuse/GROUP12 
    
fi




# Download databases
if [ ! -d ./databases ]; then
    echo "Downloading databases"
    wget_res=$(wget -q "https://www.dropbox.com/s/oqtuwkzdw5z0p6m/databases.zip")
    if [ $? -ne 0 ]
    then
        echo "Fail to download databases. The program will now terminate."
        exit 1;
    fi 
    unzip -q databases.zip
fi

mkdir -p output_defuse
Rscript --vanilla co-fuse.R 'DeFuse' './DeFuse/GROUP12' './output_defuse/recurrent' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'DeFuse' './DeFuse/GROUP1' './DeFuse/GROUP2' './output_defuse/fisher_test'

