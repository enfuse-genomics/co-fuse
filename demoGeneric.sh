 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./Generic.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://www.dropbox.com/s/m682vp6drlx8sjh/Generic.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf Generic.tar.gz
    
    #
    # Merge GROUP 1 and GROUP 2
    #
    mkdir -p Generic/GROUP12
    cp -r Generic/GROUP1/* Generic/GROUP12
    cp -r Generic/GROUP2/* Generic/GROUP12 
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

mkdir -p output_generic
Rscript --vanilla co-fuse.R 'Generic' './Generic/GROUP12' './output_generic/recurrent' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'Generic' './Generic/GROUP1' './Generic/GROUP2' './output_generic/fisher_test'


echo "Output can be found under the directory ./output_generic"

