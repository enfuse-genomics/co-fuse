 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./SoapFuse.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/SoapFuse.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf SoapFuse.tar.gz
fi

mkdir -p output_soapfuse
Rscript --vanilla co-fuse.R 'SoapFuse' './SoapFuse/GROUP1' './output_soapfuse/group1' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse.R 'SoapFuse' './SoapFuse/GROUP2' './output_soapfuse/group2' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'SoapFuse' './SoapFuse/GROUP1' './SoapFuse/GROUP2' './output_soapfuse/fisher_test'

