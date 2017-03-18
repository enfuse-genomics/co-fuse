 #!/bin/sh

TSNE_PERPLEXITY=6

if [ ! -f ./FusionCatcher.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://www.dropbox.com/s/ilfasbsr95m74a3/FusionCatcher.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf FusionCatcher.tar.gz
    
    #
    # Merge AML and MM results
    #
    mkdir -p FusionCatcher/AML_MM
    cp -r FusionCatcher/AML/* FusionCatcher/AML_MM
    cp -r FusionCatcher/MM/* FusionCatcher/AML_MM    
fi



if [ ! -d ./FusionCatcher ] && [ -f ./FusionCatcher.tar.gz ] ; then
    echo "Unzipping FusionCatcher.tar.gz"
    
    tar xzf FusionCatcher.tar.gz
    
    #
    # Merge AML and MM results
    #
    mkdir -p FusionCatcher/AML_MM
    cp -r FusionCatcher/AML/* FusionCatcher/AML_MM
    cp -r FusionCatcher/MM/* FusionCatcher/AML_MM    
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



mkdir -p output_fusioncatcher
Rscript --vanilla co-fuse.R 'FusionCatcher' './FusionCatcher/AML_MM' './output_fusioncatcher/recurrent' $TSNE_PERPLEXITY
Rscript --vanilla co-fuse2.R 'FusionCatcher' './FusionCatcher/AML' './FusionCatcher/MM' './output_fusioncatcher/fisher_test'

