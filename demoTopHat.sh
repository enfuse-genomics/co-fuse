 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./TopHat.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://www.dropbox.com/s/o8ulmwk0znzskxq/TopHat.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf TopHat.tar.gz
    
    #
    # Merge GROUP 1 and GROUP 2
    #
    mkdir -p TopHat/GROUP12
    cp -r TopHat/GROUP1/* TopHat/GROUP12
    cp -r TopHat/GROUP2/* TopHat/GROUP12     
    
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


mkdir -p output_tophat
Rscript --vanilla co-fuse.R 'TopHat' './TopHat/GROUP12' './output_tophat/recurrent' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'TopHat' './TopHat/GROUP1' './TopHat/GROUP2' './output_tophat/fisher_test'

echo "Output can be found under the directory ./output_tophat"

