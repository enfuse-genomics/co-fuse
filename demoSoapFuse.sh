 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./SoapFuse.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://www.dropbox.com/s/42eoqny3dv40zzi/SoapFuse.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf SoapFuse.tar.gz
    
    #
    # Merge GROUP 1 and GROUP 2
    #
    mkdir -p SoapFuse/GROUP12
    cp -r SoapFuse/GROUP1/* SoapFuse/GROUP12
    cp -r SoapFuse/GROUP2/* SoapFuse/GROUP12 
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


mkdir -p output_soapfuse
Rscript --vanilla co-fuse.R 'SoapFuse' './SoapFuse/GROUP12' './output_soapfuse/recurrent' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'SoapFuse' './SoapFuse/GROUP1' './SoapFuse/GROUP2' './output_soapfuse/fisher_test'

echo "Output can be found under the directory ./output_soapfuse"



